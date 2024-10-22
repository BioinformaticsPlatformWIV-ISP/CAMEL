import os
import re
import pandas as pd
from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Pharokka(Tool):
    """
    Pharokka is a fast phage genome annotation tool.
    """

    def __init__(self, camel: Camel.get_instance()) -> None:
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('Pharokka', '1.7.3', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()

        # Collect the output
        for key, basename in zip(('GBK', 'CARD', 'VFDB', 'STATS', 'INPHARED'),
                                 ('pharokka.gbk', 'top_hits_card.tsv', 'top_hits_vfdb.tsv',
                                  'pharokka_cds_functions.tsv', 'pharokka_top_hits_mash_inphared.tsv')):
            path_out = self.folder / 'pharokka' / basename
            if not path_out.exists():
                raise ToolExecutionError(f'{path_out} not generated ({key})')
            self._tool_outputs[key] = [ToolIOFile(path_out)]

        # Parse output files
        self._parse_card_hits(self._tool_outputs['CARD'][0].path)
        self._parse_vfdb_hits(self._tool_outputs['VFDB'][0].path)
        self._parse_metrics(self._tool_outputs['STATS'][0].path)

    def _check_input(self) -> None:
        """
        Checks the tool input.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise ValueError("No FASTA input found")
        if len(self._tool_inputs['FASTA']) != 1:
            raise ValueError("Only one FASTA file can be annotated at a time.")
        super(Pharokka, self)._check_input()

    def _check_command_output(self) -> None:
        """
        Checks the command output to see if the program ran correctly.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Error executing Pharokka: {}".format(self.stderr))

    def __build_command(self) -> None:
        """
        Build the command line call.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(self._build_options()),
            f'-i',
            str(self._tool_inputs['FASTA'][0].path)
        ])

    def _parse_metrics(self, tsv_path: Path) -> None:
        """
        Parses the file with the summary metrics and stores them in informs.
        :param path_csv: Path to output file
        :return: None
        """
        data = pd.read_table(tsv_path)
        cds_count = data.iloc[0]['Count']
        self._informs['stats'] = f'{cds_count} CDS were annotated by Pharokka.'

    def _parse_card_hits(self, tsv_path: Path) -> None:
        """
        Parses the top CARD hits TSV file and stores the results in the informs.
        :param path_csv: Path to output file
        :return: None
        """
        data = pd.read_table(tsv_path)
        if data.empty:
            self._informs['card_hits'] = 'No AMR gene was detected with the CARD database.'
        else:
            self._informs['card_hits'] = f'{data.shape[0]} AMR genes were detected with the CARD database.'

    def _parse_vfdb_hits(self, tsv_path: Path) -> None:
        """
        Parses the top VFDB hits TSV file and stores the results in the informs.
        :param path_csv: Path to output file
        :return: None
        """
        data = pd.read_table(tsv_path)
        if data.empty:
            self._informs['vfdb_hits'] = 'No virulence gene was detected with the VFDB database.'
        else:
            self._informs['vfdb_hits'] = f'{data.shape[0]} virulence genes were detected with the VFDB database.'
