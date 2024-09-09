from pathlib import Path
import argparse
from typing import Optional, Sequence

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool

class Bakta(Tool):
    """
    Bakta is a tool for the rapid & standardized annotation of bacterial genomes and plasmids from both isolates and MAGs
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: Camel instance
        return: None
        """
        super().__init__('Bakta', '1.9.4', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Build the Bakta command and execute it
        self.__build_command()
        self._execute_command()
        # Set the output
        self.__set_output()

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        super()._check_input()
        if 'FASTA' in self._tool_inputs:
            fasta_file = self._tool_inputs['FASTA'][0]
            self._name = fasta_file.path.stem
            if len(self._tool_inputs['FASTA']) != 1:
                raise InvalidInputSpecificationError("FASTA input requires exactly 1 file.")
        else:
            raise ValueError("FASTA input is required")

    def __build_command(self) -> None:
        """
        Builds the command.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            str(self._tool_inputs['FASTA'][0].path),
            *self._build_options()]
        )

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if 'error' in self.stderr.lower():
            raise ToolExecutionError(f"Command execution failed (stderr: {self.stderr}).")
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    def __set_output(self) -> None:
        """
        Collects the output files of interest.:
        - faa file
        - gff3 file
        - gbff file
        :return: None
        """
        self._tool_outputs['FAA_FILE'] = [ToolIOFile(self.folder / f'{self._name}.faa')]
        self._tool_outputs['GFF3_FILE'] = [ToolIOFile(self.folder / f'{self._name}.gff3')]
        self._tool_outputs['GBFF_FILE'] = [ToolIOFile(self.folder / f'{self._name}.gbff')]

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        parser.add_argument('-db', type=Path, help="Path to the database, default: /db/bakta/latest", required=True)
        parser.add_argument('-o', type=Path, help="Output directory", required=False)
        parser.add_argument('--genus', type=str, help="Genus name (i.e. Escherichia)", required=False)
        parser.add_argument('--species', type=str, help="Species name (i.e. 'coli O157:H7')", required=False)
        parser.add_argument('--proteins', type=Path, help="Fasta file of trusted protein sequences for annotation")
        parser.add_argument('--gram', type=str, help="Gram type for signal peptide predictions {+/-/?} (default = ?)", default="?")
        parser.add_argument('--compliant', action="store_true", help="Force Genbank/ENA/DDJB compliance")
        parser.add_argument('--meta', action="store_true", help="Run in metagenome mode. This only affects CDS prediction.")
        parser.add_argument('--threads', type=int, default=10, help="Number of threads")

        return parser.parse_args(args)


if __name__ == '__main__':
    main = Bakta(Camel.get_instance())


