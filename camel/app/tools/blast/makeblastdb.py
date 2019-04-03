import os

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class MakeBlastDb(Tool):
    """
    This tool can be used to create BLAST databases.
    It indexes a FASTA file in place.
    """

    def __init__(self, camel: Camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('makeblastdb', '2.6.0', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        output_file = os.path.join(self._folder, os.path.basename(self._tool_inputs['FASTA'][0].path))
        if not os.path.exists(output_file):
            os.symlink(self._tool_inputs['FASTA'][0].path, output_file)
        self.__build_command(output_file)
        self._execute_command()
        self._tool_outputs['FASTA'] = [ToolIOFile(output_file)]

    def __build_command(self, fasta_path: str) -> None:
        """
        Builds the command line call.
        :param fasta_path: FASTA path
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f'-in {fasta_path}',
        ] + self._build_options())
