import logging
import os
from collections import namedtuple

from app.components.filesystemhelper import FileSystemHelper
from app.loggers.logmanager import LogManager
from app.services.stepservice import StepService


class Step(object):
    """
    This class represents a step in a pipeline. It executes a single tool.
    """
    StepInput = namedtuple('StepInput', ('name', 'source', 'alias', 'required',))
    ExternalInput = namedtuple('ExternalInput', ('files', 'alias',))
    StepInform = namedtuple('StepInform', ('source', 'alias'))

    ConditionalNextStep = namedtuple('ConditionalNextStep', ('condition', 'step_id'))
    DefaultNextStep = namedtuple('DefaultNextStep', ('step_id',))

    def __init__(self, name, tool, pipeline, camel):
        """
        Initializes a step.
        """
        self._name = name
        self._tool = tool
        self._step_id = None
        self._input_specification = None
        self._inform_specification = None
        self._next_step_specification = None
        self._folder = None
        self._step_inputs = {}
        self._input_informs = {}
        self._pipeline = pipeline
        self._camel = camel
        self._pipeline_options = {}
        self._job_options = {}
        self._step_service = None

    def run_step(self):
        """
        Runs the current step.
        :return: None
        """
        self._create_folder()
        LogManager.attach_step_handlers(self._folder)
        self._tool.add_input_files(self._step_inputs)
        self._tool.add_input_informs(self._input_informs)
        logging.info("Default parameters loaded: {}".format(self._tool.parameter_overview))
        self._add_pipeline_parameters()
        self._add_job_parameters()
        self._tool.run(self._folder)
        logging.info('Step output: {}'.format(self.outputs.items()))
        logging.info('Step informs: {}'.format(self.informs.items()))
        if self._pipeline.job_id is not None:
            self._log_outputs()
            self._log_job_parameters()
        LogManager.detach_step_handlers()

    @property
    def name(self):
        """
        Returns the name of this step.
        :return: Name
        """
        return self._name

    @property
    def step_id(self):
        """
        Returns the database id of this step.
        :return: Database id
        """
        return self._step_id

    @step_id.setter
    def step_id(self, id_):
        """
        Sets the database id.
        :return: None
        """
        self._step_id = id_
        self._step_service = StepService(id_, self._camel.connection)

    @property
    def input_specification(self):
        """
        Returns the input specification.
        :return: Input specification
        """
        specification = []
        if self._input_specification is None:
            return specification
        for item in self._input_specification:
            try:
                if 'external-input' in item:
                    files = [item['type'](element) for element in item['external-input']]
                    specification.append(Step.ExternalInput(files, item['alias']))
                else:
                    alias = item['alias'] if 'alias' in item else item['name']
                    required = item.get('required', False)
                    specification.append(Step.StepInput(item['name'], item['from'], alias, required))
            except KeyError as err:
                raise ValueError("'{}' missing from input specification: {}".format(err, item.items()))
        return specification

    @input_specification.setter
    def input_specification(self, specification):
        """
        Sets the input specification for this step.
        :param specification: Input specification
        :return: None
        """
        self._input_specification = specification

    @property
    def inform_specification(self):
        """
        Returns the inform specification.
        :return: Inform specification
        """
        specification = []
        if self._inform_specification is None:
            return specification
        for item in self._inform_specification:
            try:
                specification.append(Step.StepInform(item['from'], item['alias']))
            except KeyError as err:
                raise ValueError("'{}' missing from inform specification: {}".format(err, item.items()))
        return specification

    @inform_specification.setter
    def inform_specification(self, specification):
        """
        Sets the inform specification for this step.
        :param specification: Inform specification
        :return: None
        """
        self._inform_specification = specification

    @property
    def next_step_specification(self):
        """
        Returns the specification of the next step.
        :return: Specification
        """
        if not self._next_step_specification:
            return None
        specification = []
        for branch in self._next_step_specification:
            try:
                if 'when' in branch:
                    specification.append(Step.ConditionalNextStep(branch['when']['condition'],
                                                                  branch['when']['step-id']))
                elif 'default' in branch:
                    specification.append(Step.DefaultNextStep(branch['default']['step-id']))
                else:
                    raise ValueError("Next step specification is not valid: {}".format(branch))
            except KeyError as err:
                raise ValueError("{} missing from step specification".format(err))
        return specification

    @next_step_specification.setter
    def next_step_specification(self, specification):
        """
        Sets the specification for the next step.
        :param specification: Specification
        :return: None
        """
        self._next_step_specification = specification

    def add_input(self, key, new_input):
        """
        Adds input to the step.
        :param key: Key for the input
        :param new_input: List of the new inputs
        :return: None
        """
        if len(new_input) == 0:
            return
        if key not in self._step_inputs:
            self._step_inputs[key] = []
        self._step_inputs[key].extend(new_input)
        logging.info("Input '{}' added: {}".format(key, new_input))

    def add_inform(self, key, inform):
        """
        Adds inform to the step.
        :param key: Key for the inform
        :param inform: Inform
        :return: None
        """
        self._input_informs[key] = inform
        logging.info("Inform '{}' added: {}".format(key, inform))

    @property
    def outputs(self):
        """
        Returns the outputs of this step.
        :return: Outputs
        """
        return self._tool.tool_outputs

    @property
    def informs(self):
        """
        Returns the informs of this step.
        :return: Informs
        """
        return self._tool.informs

    def add_pipeline_options(self, options):
        """
        Adds pipeline options to the step, pipeline options can override tool default options.
        :param options: Pipeline options
        :return: None
        """
        self._pipeline_options.update(options)

    def add_job_options(self, options):
        """
        Adds job options to the step, job options can override pipeline options.
        :param options: Job options
        :return: None
        """
        self._job_options.update(options)

    def _create_folder(self):
        """
        Creates the folder for this step.
        :return: None
        """
        self._folder = os.path.join(self._pipeline.folder, FileSystemHelper.make_valid(self.name))
        if not os.path.isdir(self._folder):
            os.makedirs(self._folder)

    def _add_pipeline_parameters(self):
        """
        Adds the pipeline parameters.
        :return: None
        """
        if len(self._pipeline_options) > 0:
            logging.info("Adding pipeline parameters")
            self._tool.update_parameters(**self._pipeline_options)

    def _add_job_parameters(self):
        """
        Adds the job parameters.
        :return: None
        """
        if len(self._job_options) > 0:
            logging.info("Adding job parameters")
            self._tool.update_parameters(**self._job_options)

    def _log_outputs(self):
        """
        Logs the outputs in the database.
        :return: None
        """
        for key, files in self.outputs.items():
            for i in range(0, len(files)):
                if files[i].logged:
                    output_data = [self._pipeline.job_id, self.step_id, files[i].TYPE_NAME, key, i, files[i].hash]
                    self._step_service.log_output(output_data)
                    logging.debug('Output {} ({}) logged'.format(key, i))

    def _log_job_parameters(self):
        """
        Logs the job parameters in the database.
        :return: None
        """
        for parameter_name, parameter_value in self._job_options.items():
            parameter_id = self._step_service.get_parameter_id(self._tool.tool_id, parameter_name)
            self._step_service.log_job_parameter(parameter_id, self.step_id, self._pipeline.job_id, parameter_value)
            logging.debug("Job parameter '{}' logged".format(parameter_name))
