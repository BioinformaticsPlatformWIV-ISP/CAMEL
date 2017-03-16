import abc
import logging

from app.command.command import Command
from app.error.invalidparametererror import InvalidParameterError
from app.error.toolexecutionerror import ToolExecutionError
from app.io.toolio import ToolIO
from app.services.toolservice import ToolService


class Tool(object):
    """
    Contains the common functionality of tools.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, name, version, camel):
        """
        Initializes a tool.
        """
        logging.debug("Initializing tool: {} {}".format(name, version))
        self._name = name
        self._version = version
        self._tool_inputs = {}
        self._tool_outputs = {}
        self._informs = {}
        self._input_informs = {}
        self._tool_service = ToolService(name, version, camel.connection)
        self._tool_command = self._tool_service.get_tool_command()
        self._parameters = self._tool_service.get_default_parameters()
        self._command = Command()
        self._folder = None

    @property
    def name(self):
        """
        Returns the name of this tool.
        :return: Name
        """
        return '{} {}'.format(self._name, self._version)

    @property
    def tool_id(self):
        """
        Returns the tool id.
        :return: Tool id
        """
        return self._tool_service.tool_id

    @property
    def tool_outputs(self):
        """
        Returns the tool outputs.
        :return: Tool outputs
        """
        return self._tool_outputs

    @property
    def informs(self):
        """
        Returns the tool informs.
        :return: Informs
        """
        return self._informs

    @property
    def stdout(self):
        """
        Returns the command line stdout.
        :return: Stdout
        """
        return self._command.stdout

    @property
    def stderr(self):
        """
        Returns the command line stderr.
        :return: Stderr
        """
        return self._command.stderr

    @property
    def parameter_overview(self):
        """
        Returns an overview of the parameters as a string.
        :return: Parameters overview
        """
        return ', '.join(["{}: '{}'".format(p, self._parameters[p].value) for p in sorted(self._parameters)]) if \
            len(self._parameters) > 0 else '/'

    def add_input_files(self, input_files):
        """
        Updates the input files for a tool.
        :param input_files: New input files
        :return: None
        """
        for key, items in input_files.iteritems():
            if key in self._tool_inputs:
                self._tool_inputs[key] += items
            else:
                self._tool_inputs[key] = items

    def add_input_informs(self, informs):
        """
        Updates the input informs for a tool.
        :param informs: New informs
        :return: None
        """
        self._input_informs.update(informs)

    def update_parameters(self, **kwargs):
        """
        Updates the parameters for this tool.
        :param kwargs: Parameters in key value format
        :return: None
        """
        for parameter_name, new_value in kwargs.iteritems():
            parameter = self._tool_service.get_parameter(parameter_name)
            if not parameter:
                raise InvalidParameterError("{} has no parameter '{}'".format(self._name, parameter_name))
            if new_value is False:
                if parameter_name not in self._parameters:
                    raise ValueError("Cannot disable parameter '{}' (not present in parameters)".format(parameter_name))
                logging.info("Disabling parameter: {}".format(parameter_name))
                del(self._parameters[parameter_name])
            else:
                if new_value is True or new_value is None:
                    parameter.value = None
                else:
                    parameter.value = str(new_value)
                if parameter_name not in self._parameters:
                    logging.info("Parameter '{}' added, value: {}".format(parameter_name, parameter.value))
                else:
                    old_value = self._parameters[parameter_name].value
                    logging.info("Parameter '{}' value '{}' changed to '{}'".format(
                        parameter_name, old_value, new_value))
                self._parameters[parameter_name] = parameter

    def clear_parameters(self):
        """
        Clears all the parameters of the given tool.
        :return: None
        """
        logging.info("Removing {} parameters".format(len(self._parameters)))
        self._parameters.clear()

    def run(self, folder='.'):
        """
        Runs this tool.
        :param folder: Folder to run the tool in.
        :return: None
        """
        self._folder = folder
        logging.info('Running tool {}'.format(self.name))
        logging.info('Working directory: {}'.format(self._folder))
        logging.info('Tool parameters: {}'.format(self.parameter_overview))
        self._check_parameters()
        self._check_input()
        self._execute_tool()
        self._check_output()

    def get_outputs(self, key):
        """
        Returns the outputs with the given key.
        :param key: output key
        :return: Output list
        """
        if key not in self._tool_outputs:
            raise ValueError("No output file with key '{}' found".format(key))
        return self._tool_outputs[key]

    def _build_dependencies(self):
        """
        Builds the dependencies.
        :return: Command to load dependencies
        """
        dependencies = self._tool_service.get_dependencies()
        if dependencies is not None:
            return 'module load ' + dependencies + '; '
        else:
            return ''

    def _build_options(self, excluded_parameters=None, delimiter=' '):
        """
        Builds the options string.
        :parameter delimiter: Delimiter between option and value
        :return: Options string
        """
        options = []
        for name, parameter in self._parameters.iteritems():
            if (excluded_parameters is not None) and (name in excluded_parameters):
                continue
            if parameter.value is not None:
                options.append(parameter.option + delimiter + str(parameter.value))
            else:
                options.append(parameter.option)
        return options

    def _execute_command(self, folder=None):
        """
        Executes a the command.
        :return: None
        """
        if folder is None:
            folder = self._folder
        if self._command.command is None:
            raise ValueError("Command is 'None'.")
        self._command.command = self._build_dependencies() + self._command.command
        self._command.run_command(folder)
        self._check_command_output()

    def _check_command_output(self):
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if self.stderr != '':
            raise ToolExecutionError("Command execution failed (stderr: {}).".format(self.stderr))
        elif self._command.returncode != 0:
            raise ToolExecutionError("Command execution failed (Exit code: {})".format(self._command.returncode))

    @abc.abstractmethod
    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        pass

    def _check_parameters(self):
        """
        Checks if the tool parameters are valid.
        :return: None
        """
        mandatory_parameters = self._tool_service.get_names_mandatory_parameter()
        for mandatory_parameter in mandatory_parameters:
            if mandatory_parameter not in self._parameters:
                raise ValueError("Mandatory parameter {} not set".format(mandatory_parameter))

    def _check_input(self):
        """
        Checks if the tool input is valid.
        :return: None
        """
        for input_key, input_list in self._tool_inputs.iteritems():
            for tool_input in input_list:
                if tool_input is None:
                    raise ValueError("Tool input with key {} is None".format(input_key))
                if not tool_input.is_valid():
                    raise ValueError("Invalid tool input with key {}: {}".format(input_key, tool_input))

    def _check_output(self):
        """
        Checks if the output is valid.
        :return: None
        """
        for output_key, output_list in self._tool_outputs.iteritems():
            for tool_output in output_list:
                if tool_output is None:
                    raise ValueError("Tool output with key {} is None".format(output_key))
                if not isinstance(tool_output, ToolIO):
                    raise ValueError("'{} {}' is not a tool output object".format(tool_output, type(tool_output)))
                if not tool_output.is_valid():
                    raise ValueError("Invalid tool output with key {}: {}".format(output_key, tool_output))
