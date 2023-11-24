from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftoolsbase import BcftoolsBase


class BcftoolsConsensus(BcftoolsBase):
    """
    Create consensus sequences by applying VCF variants to a reference fasta file.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('bcftools consensus', '1.17', camel)

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if not any(x in self._tool_inputs for x in ('VCF', 'VCF_GZ')):
            raise InvalidInputSpecificationError("No variants input found (VCF/VCF_GZ).")
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No reference FASTA input found.")
        super(BcftoolsConsensus, self)._check_input()

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
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f'--fasta-ref {self._tool_inputs["FASTA"][0].path}',
            str(next(self._tool_inputs[k][0].path for k in ('VCF', 'VCF_GZ') if k in self._tool_inputs)),
            ' '.join(self._build_options())
        ])

    def __set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['FASTA'] = [ToolIOFile(self._get_output_path())]
