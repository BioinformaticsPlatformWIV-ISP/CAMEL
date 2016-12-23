import logging
import os

import time
import yaml

from app.components.filesystemhelper import FileSystemHelper
from step import Step


class Pipeline(object):
    """
    Class meant to handle the workflow of steps.
    """

    def __init__(self, yaml_file, camel):
        """
        Initializes a pipeline.
        :param yaml_file: Pipeline YAML file
        :param camel: CAMEL instance
        """
        self._camel = camel
        self._name = None
        self._initial_input = None
        self._steps = None
        self._destination_path = None
        self._folder = None
        self._parse_yaml_file(yaml_file)

    @property
    def folder(self):
        """
        Returns the output folder of this pipeline.
        :return: Output folder
        """
        return self._folder

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
        stack = []
        stack.insert(0, self._steps[0])
        self._steps[0].step_inputs = self._initial_input
        while stack:
            logging.info('Running step {}'.format(stack[0].name))
            self._prepare_inputs_for_step(stack[0])
            stack[0].run_step()
            next_step = self._get_next_step(stack[0])
            if next_step is None:
                break
            stack.insert(1, next_step)
            stack.pop(0)
        logging.info('Finished running pipeline')

    def add_options(self, options):
        """
        Adds options that override the default tool options.
        :param options: Options with step name as key
        :return: None
        """
        for step_name, step_options in options.iteritems():
            step = self._get_step(step_name)
            if step is None:
                raise ValueError("No step named '{}'".format(step_name))
            step.add_options(step_options)

    def _parse_yaml_file(self, yaml_file):
        """
        Parses the YAML file to get the name and the steps.
        :param yaml_file: YAML file
        :return: None
        """
        with open(yaml_file) as input_handle:
            yaml_data = yaml.load(input_handle)
        try:
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
        self._steps = []
        for step in step_data:
            try:
                new_step = Step(step['id'], step['tool'](self._camel), step['next-step'], self)
                if 'inputs' in step:
                    new_step.input_specification = step['inputs']
                self._steps.append(new_step)
            except KeyError as err:
                raise ValueError("'{}' missing in step specification {}".format(err.message, step_data.index(step)))

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
        # Todo default always next step in yaml?
        next_step_name = current_step._next_step[0]['default']['step-id']
        if next_step_name == 'exit':
            return None
        next_step = self._get_step(next_step_name)
        if next_step is None:
            raise ValueError("Step '{}' is missing in step specification".format(next_step_name))
        return next_step

    def _prepare_inputs_for_step(self, step):
        """
        Returns the inputs for the given step.
        :param step: Step
        :return: None
        """
        for input_ in step.input_specification:
            logging.info('Preparing input: {}'.format(str(input_)))
            if type(input_) is Step.StepInput:
                if input_.source == 'initial':
                    step.add_input(input_.alias, self._initial_input[input_.name])
                elif self._get_step(input_.source) is not None:
                    step.add_input(input_.alias, self._get_step(input_.source).outputs[input_.name])
                else:
                    raise ValueError("No step named '{}'".format(input_.source))
            if type(input_) is Step.ExternalInput:
                step.add_input(input_.alias, input_.files)

    def _get_step(self, step_name):
        """
        Returns the step with the given name.
        :param step_name: Step name
        :return: Step
        """
        for step in self._steps:
            if step.name == step_name:
                return step
        return None
