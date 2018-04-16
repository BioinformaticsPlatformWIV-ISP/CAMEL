import logging
from typing import Dict, List

from snakemake.io import Wildcards

from camel.app.camel import Camel
from camel.app.io.toolio import ToolIO
from camel.app.loggers.logmanager import LogManager
from camel.app.services.stepservice import StepService
from camel.app.tools.tool import Tool


class Step(object):
    """
    This class represents a step in a Snakemake pipeline. It executes a single tool.
    """

    def __init__(self, rule_name: str, tool: Tool, camel: Camel, folder: str, config: dict,
                 wildcards: Wildcards=None, pipeline_output: bool=False, log_step: bool=None) -> None:
        """
        Initializes a step.
        :param rule_name: Name of the snakerule
        :param tool: Tool object of the tool that needs to be executed
        :param camel: Camel object
        :param folder: Folder in which the step is being run
        :param config: Snakemake config dictionary
        :param wildcards: Wildcards object from snakemake
        :param pipeline_output: Boolean to indicate whether outputs are pipeline outputs
        :param log_step: Boolean to indicate whether outputs for this step have to be logged
        """
        self._name = rule_name
        self._tool = tool
        self._step_inputs = {}
        self._input_informs = {}
        self._camel = camel
        self._pipeline_options = {}
        self._job_options = {}
        self._db_logging = Step.step_is_logged(config.get('logging_level'), log_step, pipeline_output)
        self._step_service = StepService(camel.connection)
        self._folder = folder
        self._job_id = config['pipeline_job_id'] if self._db_logging else None
        self._pipeline_output = pipeline_output
        self._wildcards = wildcards

    @property
    def name(self) -> str:
        """
        Returns the name of this step.
        :return: Name
        """
        return self._name

    @property
    def outputs(self) -> Dict[str, List[ToolIO]]:
        """
        Returns the outputs of this step.
        :return: Outputs
        """
        return self._tool.tool_outputs

    @property
    def informs(self) -> dict:
        """
        Returns the informs of this step.
        :return: Informs
        """
        return self._tool.informs

    def add_inputs(self, dict_: dict) -> None:
        """
        Adds the inputs to the step
        :param dict_: Dictionary with input objects
        :return: None
        """
        self._step_inputs = dict_

    def add_informs(self, dict_: dict) -> None:
        """
        Adds informs to the step.
        :param dict_: Dictionary with the informs
        :return: None
        """
        self._input_informs = dict_
        logging.info("Inform added: {}".format(dict_))

    def add_job_options(self, options: dict) -> None:
        """
        Adds job options to the step, job options can override pipeline options.
        :param options: Job options
        :return: None
        """
        self._job_options.update(options)

    def _add_job_parameters(self) -> None:
        """
        Adds the job parameters.
        :return: None
        """
        if len(self._job_options) > 0:
            logging.info("Adding job parameters")
            self._tool.update_parameters(**self._job_options)

    def run_step(self) -> None:
        """
        Runs the current step.
        :return: None
        """
        LogManager.attach_step_handlers(self._folder)
        self._tool.add_input_files(self._step_inputs)
        self._tool.add_input_informs(self._input_informs)
        logging.info("Default parameters loaded: {}".format(self._tool.parameter_overview))
        self._add_job_parameters()
        self._tool.run(self._folder)
        logging.info('Step output: {}'.format(list(self.outputs.items())))
        logging.info('Step informs: {}'.format(list(self.informs.items())))
        if self._db_logging:
            self._log_outputs()
        LogManager.detach_step_handlers()

    def _log_outputs(self) -> None:
        """
        Logs the outputs in the database.
        :return: None
        """
        logging.info("Logging step outputs")
        for key, files in self.outputs.items():
            for i in range(0, len(files)):
                if files[i].logged:
                    output_data = (self._job_id, self._name, self._wildcards, files[i].type_name, key, i,
                                   files[i].hash, self._pipeline_output)
                    logging.debug('OUTPUT DATA: {}'.format(output_data))
                    self._step_service.log_output(output_data)
                    logging.debug('Output {} ({}) logged'.format(key, i))

    @staticmethod
    def step_is_logged(logging_level: str, log_step: bool, pipeline_output: bool) -> bool:
        """
        This helper function is used to check whether a step should be logged depending on the logging level set in the
        config and the variables passed to the step object.
        :param logging_level: Logging level as defined in the config
        :param log_step: Boolean to indicate if this step should be logged (overwrites the logging level is it is not
        :param pipeline_output: True if this step is a pipeline output
        set to None.
        :return: True if the step should be logged, False otherwise.
        """
        if logging_level == 'step':
            if log_step is None or log_step is True:
                return True
            else:
                return False
        elif logging_level == 'pipeline':
            if log_step is True or pipeline_output is True:
                return True
            else:
                return False
        else:
            return False
