import logging

from app.pipeline.step import Step
from app.services.pipelineservice import PipelineService
from app.services.stepservice import StepService


class SnakeStep(Step):
    """
    This class represents a step in a Snakemake pipeline. It executes a single tool.
    """

    def __init__(self, name, tool, camel, folder, config):
        """
        Initializes the step
        :param name: Name of the step
        :param tool: Tool that needs to be executed
        :param config: Dictionary with a Camel instance, logging (boolean), pipeline_name and job_id values 
        :param folder: Working directory where the tool has to run
        """
        if not all([x in config for x in ['logging', 'pipeline_name', 'pipeline_job_id']]):
            raise ValueError('Not all required keys (logging, pipeline_name, pipeline_job_id) are given in the config dictionary: {}'.format(config))
        super(SnakeStep, self).__init__(name, tool, None, camel)
        self._db_logging = config['logging']
        self._pipeline_service = PipelineService(config['pipeline_name'], camel.connection)
        self.__set_step_id()
        self._folder = folder
        self._job_id = config['pipeline_job_id']

    def __set_step_id(self):
        """
        Sets the step id when logging is requested in the pipeline (= a pipeline service is created)
        :return: None
        """
        if self._db_logging is True:
            self.step_id = self._pipeline_service.get_step_id(self.name)
            if self.step_id is None:
                self.step_id = self._pipeline_service.insert_step(self.name)

    def run_step(self):
        """
        Runs the current step.
        :return: None
        """
        self._tool.add_input_files(self._step_inputs)
        self._tool.add_input_informs(self._input_informs)
        logging.info("Default parameters loaded: {}".format(self._tool.parameter_overview))
        self.__set_pipeline_options()
        self._add_pipeline_parameters()
        self._add_job_parameters()
        self._tool.run(self._folder)
        logging.info('Step output: {}'.format(list(self.outputs.items())))
        logging.info('Step informs: {}'.format(list(self.informs.items())))
        if self._db_logging is True:
            self._log_outputs()
            self._log_job_parameters()

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

    def add_inputs(self, dict_):
        """
        Adds the inputs to the step
        :param dict_: Dictionary with input objects
        :return: None
        """
        self._step_inputs = dict_

    def add_informs(self, dict_):
        """
        Adds informs to the step.
        :param dict_: Dictionary with the informs
        :return: None
        """
        self._input_informs = dict_
        logging.info("Inform added: {}".format(dict_))

    def __set_pipeline_options(self):
        """
        Loads options from the database.
        :return: None
        """
        logging.info("Loading pipeline parameters from database")
        if self._pipeline_service is not None:
            for parameter in self._pipeline_service.get_pipeline_step_parameters(self.name):
                self.add_pipeline_options({parameter.name: parameter.value})

    def _log_outputs(self):
        """
        Logs the outputs in the database.
        :return: None
        """
        for key, files in self.outputs.items():
            for i in range(0, len(files)):
                if files[i].logged:
                    output_data = [self._job_id, self.step_id, files[i].TYPE_NAME, key, i, files[i].hash]
                    logging.debug('OUTPUT DATA: {}'.format(output_data))
                    self._step_service.log_output(output_data)
                    logging.debug('Output {} ({}) logged'.format(key, i))

    def _log_job_parameters(self):
        """
        Logs the job parameters in the database.
        :return: None
        """
        for parameter_name, parameter_value in self._job_options.items():
            parameter_id = self._step_service.get_parameter_id(self._tool.tool_id, parameter_name)
            self._step_service.log_job_parameter(parameter_id, self.step_id, self._job_id, parameter_value)
            logging.debug("Job parameter '{}' logged".format(parameter_name))
