import logging

from abc import ABC, abstractmethod

from app.services.stepservice import StepService


class Step(ABC):
    """
    This class represents a step in a pipeline. It executes a single tool.
    """

    def __init__(self, name, tool, pipeline, camel):
        """
        Initializes a step.
        :param name: Name of the step
        :param tool: Tool that needs to be executed
        :param pipeline: Pipeline object of which this step is a part
        :param camel: Camel object
        """
        self._name = name
        self._tool = tool
        self._step_id = None
        self._folder = None
        self._step_inputs = {}
        self._input_informs = {}
        self._pipeline = pipeline
        self._camel = camel
        self._pipeline_options = {}
        self._job_options = {}
        self._step_service = None

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

    @abstractmethod
    def run_step(self):
        pass
