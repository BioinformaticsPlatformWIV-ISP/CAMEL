from pathlib import Path
from statistics import stdev

from Bio import SeqIO

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Prodigal(Tool):
    """
    Fast, reliable protein-coding gene prediction for prokaryotic genomes.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('Prodigal', '2.6.3', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        output_gbk = self.folder / 'coords.gbk'
        output_fasta = self.folder / 'protein_sequences.faa'
        self._command.command = ' '.join([
            self._tool_command,
            f"-i {self._tool_inputs['FASTA'][0].path}",
            f'-o {output_gbk}',
            f'-a {output_fasta}',
            *self._build_options(),
        ])
        self._execute_command()
        self._tool_outputs['FASTA'] = [ToolIOFile(output_fasta)]
        self._tool_outputs['GBK'] = [ToolIOFile(output_gbk)]
        self._extract_fasta_stats(output_fasta)

    def _check_input(self) -> None:
        """
        Checks the tool input.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("FASTA input is required")
        super()._check_input()

    def _check_command_output(self) -> None:
        """
        Checks the command output to see if the program ran correctly.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f'Error executing {self.name}: {self.stderr}')

    def _extract_fasta_stats(self, path_fasta: Path) -> None:
        """
        Extracts stats from the FASTA output file and stores them in the informs.
        :param path_fasta: Input FASTA file
        :return: None
        """
        with path_fasta.open() as handle:
            seqs = list(SeqIO.parse(handle, 'fasta'))
        self._informs['cds'] = {
            'nb': len(seqs),
            'avg_len': sum(len(s.seq) for s in seqs) / len(seqs),
            'std': stdev([len(s.seq) for s in seqs])
        }
