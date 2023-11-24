from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftoolsbase import BcftoolsBase


class BcftoolsView(BcftoolsBase):
    """
    VCF/BCF conversion, view, subset and filter VCF/BCF files.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('bcftools view', '1.17', camel)

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('BCF', 'VCF', 'VCF_GZ')):
            raise InvalidInputSpecificationError("No input file found (BCF / VCF_GZ / VCF supported)")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __build_command(self) -> None:
        """
        Builds the command.
        :return: None
        """
        command_parts = [
            self._tool_command,
            str(next(
                self._tool_inputs[k][0].path for k in ('VCF', 'VCF_GZ', 'BCF') if k in self._tool_inputs))
        ]
        command_parts += self._build_options(['compress_output'])
        self._command.command = ' '.join(command_parts)

    def __set_output(self) -> None:
        """
        Sets the tool output.
        :return: None
        """
        self._tool_outputs[self._get_output_key()] = [ToolIOFile(self._get_output_path())]
