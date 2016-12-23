import logging
import os
from collections import namedtuple

from app.components.filesystemhelper import FileSystemHelper


class Step(object):
    """
    This class represents a step in a pipeline. It executes a single tool.
    """
    StepInput = namedtuple('StepInput', ('name', 'source', 'alias',))
    ExternalInput = namedtuple('ExternalInput', ('files', 'alias',))

    def __init__(self, name, tool, next_step, pipeline):
        """
        Initializes a step.
        """
        self._name = name
        self._tool = tool
        self._next_step = next_step
        self._input_specification = None
        self._folder = None
        self._step_inputs = {}
        self._pipeline = pipeline
        self._options = {}

    @property
    def name(self):
        """
        Returns the name of this step.
        :return: Name
        """
        return self._name

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
                    specification.append(Step.StepInput(item['name'], item['from'], alias))
            except KeyError as err:
                raise ValueError("'{}' missing from input specification: {}".format(err.message, item.items()))
        return specification

    @input_specification.setter
    def input_specification(self, specification):
        """
        Sets the input specification for this step.
        :param specification: Input specification
        :return: None
        """
        self._input_specification = specification

    def add_input(self, key, new_input):
        """
        Adds input to the step.
        :param key: Key for the input
        :param new_input: List of the new inputs
        :return: None
        """
        logging.debug("Input '{}' added: {}".format(key, new_input))
        if len(new_input) == 0:
            return
        if key in self._step_inputs:
            self._step_inputs[key].extend(new_input)
        else:
            self._step_inputs[key] = new_input

    @property
    def outputs(self):
        """
        Returns the outputs of this step.
        :return: Outputs
        """
        return self._tool._tool_outputs

    def add_options(self, options):
        """
        Adds options to the step.
        :param options: Step options
        :return: None
        """
        self._options.update(options)

    def _create_folder(self):
        """
        Creates the folder for this step.
        :return: None
        """
        self._folder = os.path.join(self._pipeline.folder, FileSystemHelper.make_valid(self.name))
        if not os.path.isdir(self._folder):
            os.makedirs(self._folder)

    def run_step(self):
        """
        Runs the current step.
        :return: None
        """
        self._create_folder()
        self._tool.add_input_files(self._step_inputs)
        self._tool.update_parameters(**self._options)
        self._tool.run(self._folder)
        logging.debug('Step output: {}'.format(self.outputs.items()))
