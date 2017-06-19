import shutil
from tempfile import NamedTemporaryFile

import os

from app.components.images.svgconvert import SVGConvert
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.invalidparametererror import InvalidParameterError
from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class TreeVector(Tool):
    """
    TreeVector is a utility to create and integrate phylogenetic trees as Scalable Vector Graphics (SVG) files.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(TreeVector, self).__init__('TreeVector', '1.0', camel)

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if 'NWK' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No Newick tree input found")
        super(TreeVector, self)._check_input()

    def _check_parameters(self):
        """
        Checks if the parameters are valid.
        :return: None
        """
        if self._parameters['type'].value not in ('clad', 'simpleclad', 'phylo'):
            raise InvalidParameterError("Invalid parameter type value: {}".format(self._parameters['type'].value))
        if self._parameters['output_format'].value not in ('svg', 'png'):
            raise InvalidParameterError("Output format must be either 'png' or 'svg'")
        super(TreeVector, self)._check_parameters()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()
        if self._parameters['output_format'].value == 'png':
            self.__convert_output()

    def __build_command(self):
        """
        Builds the command line call.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            self._tool_inputs['NWK'][0].path,
            '-{}'.format(self._parameters['type'].value),
            ' '.join(self._build_options(excluded_parameters=['type', 'output_format']))
        ])

    def _check_command_output(self):
        """
        Checks if the command ran successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Error executing TreeVector: {}".format(self._command.stderr))

    def __set_output(self):
        """
        Sets the tool output.
        :return: None
        """
        self._tool_outputs['SVG'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output_filename'].value))]

    def __convert_output(self):
        """
        Converts the tool output to PNG.
        :return: None
        """
        output_file = self._tool_outputs.pop('SVG')[0].path
        temp_file = NamedTemporaryFile()
        shutil.move(output_file, temp_file.name)
        SVGConvert.convert_svg(temp_file.name, output_file)
        self._tool_outputs['PNG'] = [ToolIOFile(output_file)]
