import logging
import time

import os
import re
import yaml

from app.components.filesystemhelper import FileSystemHelper
from app.loggers.logmanager import LogManager
from app.services.pipelineservice import PipelineService
from step import Step


class Pipeline(object):
    """
    Class meant to handle the workflow of steps.
    """

    def __init__(self, yaml_files, camel, db_pipeline_parameters=False, db_logging=False):
        """
        Initializes a pipeline.
        :param yaml_files: Pipeline YAML files
        :param camel: CAMEL instance
        :param db_pipeline_parameters: Use pipeline parameters from the database.
        :param db_logging: If True, inputs & outputs are logged in the database.
        """
        self._camel = camel
        self._db_logging = db_logging
        self._name = None
        self._initial_input = None
        self._steps = []
        self._destination_path = None
        self._folder = None
        self._parse_yaml_files(yaml_files)
        self._pipeline_service = None
        self._job_id = None
        if db_pipeline_parameters is True:
            self._pipeline_service = PipelineService(self._name, camel.connection)
            self._update_steps_in_db()
            self._load_db_pipeline_parameters()
        if db_logging is True:
            if db_pipeline_parameters is False:
                raise ValueError("Cannot log outputs if 'db_parameters' is False.")
            self._pipeline_service = PipelineService(self._name, camel.connection)
            self._job_id = self._pipeline_service.insert_pipeline_job()

    @property
    def folder(self):
        """
        Returns the output folder of this pipeline.
        :return: Output folder
        """
        return self._folder

    @property
    def job_id(self):
        """
        Returns the job id of this pipeline.
        :return: Job id
        """
        return self._job_id

    @property
    def steps(self):
        """
        Returns the pipeline steps.
        :return: List of step objects
        """
        return self._steps

    def set_initial_input(self, files):
        """
        Sets the initial input files of the pipeline.
        :param files: dictionary of files to import
        :return: None
        """
        if type(files) is dict:
            self._initial_input = files
        else:
            raise TypeError("Input object should be a dictionary")

    def run(self, destination_path):
        """
        Runs the pipeline.
        :return: None
        """
        self._create_folder(destination_path)
        logging.info("Working directory: {}".format(self._folder))
        LogManager.attach_pipeline_handlers(self._folder)
        logging.info("Running pipeline {}".format(self._name))
        self._steps[0].step_inputs = self._initial_input
        if self._db_logging:
            self._log_initial_input()
        self._execute_pipeline_steps()
        logging.info("Finished running pipeline")
        LogManager.detach_pipeline_handlers()

    def _execute_pipeline_steps(self):
        """
        Executes the steps in the pipeline.
        :return: None
        """
        stack = []
        stack.insert(0, self._steps[0])
        while stack:
            logging.info('Running step {}'.format(stack[0].name))
            self._prepare_inputs_for_step(stack[0])
            self._prepare_informs_for_step(stack[0])
            stack[0].run_step()
            next_step = self._get_next_step(stack[0])
            if next_step is None:
                break
            stack.insert(1, next_step)
            stack.pop(0)

    def add_job_options(self, options):
        """
        Adds options that override the pipeline options.
        :param options: Options with step name as key
        :return: None
        """
        for step_name, step_options in options.iteritems():
            step = self.get_step(step_name)
            if step is None:
                raise ValueError("No step named '{}'".format(step_name))
            step.add_job_options(step_options)

    def _parse_yaml_files(self, yaml_files):
        """
        Parses the YAML file to get the name and the steps.
        :param yaml_files: list of YAML files
        :return: None
        """
        for yaml_file in yaml_files:
            with open(yaml_file) as input_handle:
                yaml_data = yaml.load(input_handle)
                try:
                    if self._name is None:
                        self._name = yaml_data['name']
                    self._parse_steps(yaml_data['steps'])
                except KeyError as err:
                    raise ValueError("'{}' missing in pipeline specification '{}'".format(err.message, yaml_file))

    def _parse_steps(self, step_data):
        """
        Generates steps based on the given step data structure.
        :param step_data: Step data
        :return: None
        """
        for step in step_data:
            try:
                new_step = Step(step['id'], step['tool'](self._camel), self, self._camel)
                if 'inputs' in step:
                    new_step.input_specification = step['inputs']
                if 'informs' in step:
                    new_step.inform_specification = step['informs']
                if 'next-step' in step:
                    new_step.next_step_specification = step['next-step']
                self._steps.append(new_step)
            except KeyError as err:
                raise ValueError("'{}' missing in step specification {}".format(err.message, step_data.index(step)))

    def _update_steps_in_db(self):
        """
        Updates the pipeline steps in the database. Reports a warning if a step in the database in not present in the
        YAML file.
        :return: None
        """
        step_ids = self._pipeline_service.get_step_ids()
        for step in self._steps:
            if step.name not in step_ids:
                step.step_id = self._pipeline_service.insert_step(step.name)
            else:
                step.step_id = step_ids.pop(step.name)
        if len(step_ids) > 0:
            logging.warning('Step(s) found in database that are not in YAML specification: {}'.format(
                ', '.join(step_ids.keys())))

    def _load_db_pipeline_parameters(self):
        """
        Loads options from the database.
        :return: None
        """
        logging.info("Loading pipeline parameters from database")
        for step_name, parameter in self._pipeline_service.get_pipeline_parameters():
            step = self.get_step(step_name)
            if step is None:
                raise ValueError("Cannot add parameter '{}', step '{}' does not exist".format(
                    parameter.name, step_name))
            step.add_pipeline_options({parameter.name: parameter.value})

    def _create_folder(self, destination_path):
        """
        Creates tree structure to store the pipeline results.
        pipeline name > date > time
        :param destination_path: Destination path of the pipeline
        :return: None
        """
        self._destination_path = destination_path
        date_path = time.strftime('%Y-%m-%d_%H%M%S')
        self._folder = os.path.join(destination_path, FileSystemHelper.make_valid(self._name), date_path)
        if not os.path.exists(self._folder):
            os.makedirs(self._folder)

    def _get_next_step(self, current_step):
        """
        Returns the step that needs to be executed after the current step.
        :return: Next step
        """
        next_step_name = None
        if current_step.next_step_specification is None:
            index = self._steps.index(current_step)
            try:
                next_step_name = self._steps[index+1].name
            except IndexError:
                next_step_name = 'exit'
        else:
            for branch in current_step.next_step_specification:
                if type(branch) is Step.ConditionalNextStep:
                    logging.info("Evaluating condition: {}".format(branch.condition))
                    condition_expression_resolved = self._resolve_variables(branch.condition, current_step)
                    logging.info("Condition after replacing variables: {}".format(condition_expression_resolved))
                    if self._evaluate_expression(condition_expression_resolved):
                        logging.info('Condition evaluated to True, next step: {}'.format(branch.step_id))
                        next_step_name = branch.step_id
                        break
                    else:
                        logging.info('Condition evaluated to False')
                elif type(branch) is Step.DefaultNextStep:
                    next_step_name = branch.step_id
                    break
                else:
                    raise ValueError("Cannot determine next step")

        if next_step_name == 'exit':
            return None
        next_step = self.get_step(next_step_name)
        if next_step is None:
            raise ValueError("Step '{}' is missing in step specification".format(next_step_name))
        return next_step

    def _prepare_inputs_for_step(self, step):
        """
        Prepares the inputs for the given step.
        :param step: Step
        :return: None
        """
        for input_ in step.input_specification:
            logging.info("Preparing input: {}".format(str(input_)))
            if type(input_) is Step.StepInput:
                if input_.source == 'initial':
                    try:
                        step.add_input(input_.alias, self._initial_input[input_.name])
                    except KeyError:
                        raise ValueError("Initial input does not contain key '{}'".format(input_.name))
                elif self.get_step(input_.source) is not None:
                    try:
                        step.add_input(input_.alias, self.get_step(input_.source).outputs[input_.name])
                    except KeyError as err:
                        logging.warning("Step {} has no output {}".format(input_.source, err.message))
                else:
                    raise ValueError("No step named '{}'".format(input_.source))
            elif type(input_) is Step.ExternalInput:
                step.add_input(input_.alias, input_.files)

    def _prepare_informs_for_step(self, step):
        """
        Prepares the input informs for the given step.
        :param step: Step
        :return: None
        """
        for inform in step.inform_specification:
            logging.info("Preparing inform: {}".format(str(inform)))
            source_step = self.get_step(inform.source)
            if source_step is None:
                raise ValueError("No step named '{}'".format(inform.source))
            step.add_inform(inform.alias, source_step.informs)

    def get_step(self, step_name):
        """
        Returns the step with the given name.
        :param step_name: Step name
        :return: Step
        """
        for step in self._steps:
            if step.name == step_name:
                return step
        return None

    def _resolve_variables(self, expression, current_step):
        """
        Replaces the variables in the input by their value.
        :param expression: Input expression.
        :param current_step: The current step
        :return: String with resolved variables
        """
        all_matches = re.findall(r'(\$[\w\\.]+)', expression)
        for match in all_matches:
            if match.count('.') == 1:
                step_name, key = match[1:].split('.')
                expression = expression.replace(match, '{!r}'.format(self._get_inform_value(key, step_name)))
            elif match.count('.') == 0:
                key = match[1:]
                expression = expression.replace(match, '{!r}'.format(self._get_inform_value(key, current_step.name)))
            else:
                raise ValueError("Invalid condition: {}".format(match))
        return expression

    def _get_inform_value(self, key, step_name):
        """
        Returns the given inform.
        :param key: Key
        :param step_name: Step name
        :return: Inform value
        """
        try:
            return self.get_step(step_name).informs[key]
        except KeyError as err:
            raise ValueError("Cannot retrieve '{}' informs from '{}'".format(err.message, step_name))

    @staticmethod
    def _evaluate_expression(expression):
        """
        Evaluates a string expression.
        :param expression: Expression to evaluate
        :return: Return value of the expression.
        """
        try:
            result = eval(expression)
        except SyntaxError as err:
            raise StandardError("Invalid condition: {0} ({1})".format(expression, err))

        if not type(result) is bool:
            raise StandardError('Condition did not evaluate to boolean value ({0})'.format(expression))
        return result

    def _log_initial_input(self):
        """
        Logs the initial input of the pipeline.
        :return: None
        """
        for key, files in self._initial_input.iteritems():
            for i in range(0, len(files)):
                if files[i].logged:
                    self._pipeline_service.log_initial_input(self._job_id, files[i].TYPE_NAME, key, i, files[i].hash)
                    logging.debug('Initial input {} ({}) logged'.format(key, i))
