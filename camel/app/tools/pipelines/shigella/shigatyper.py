import logging
import shutil
import pandas as pd

from pathlib import Path

from camel.app.camel import Camel
from camel.app.tools.tool import Tool
from camel.app.io.tooliofile import ToolIOFile
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError


class ShigaTyper(Tool):
    """
    ShigaTyper does Shigella/EIEC identification and serotyping based on WGS data.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the ShigaTyper tool.
        :param camel: CAMEL instance
        """
        super().__init__('ShigaTyper', '2.0.5', camel)

    def _check_input(self) -> None:
        """
        Checks whether the provided input is valid:
        - Illumina paired-end reads are the only required input
        :return: None
        """
        if 'FASTQ_FWD' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Paired-end reads are required')
        if 'FASTQ_REV' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Paired-end reads are required')
        super()._check_input()

    def __build_command(self, input_fwd: Path, input_rev: Path, sample_name: str) -> None:
        """
        Concatenates required parameters and options to build the command.
        :param input_fwd: Path to forward fastq
        :param input_rev: Path to reverse fastq
        :param sample_name: Sample ID
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f'--R1 {input_fwd}',
            f'--R2 {input_rev}',
            f'--name {sample_name}',
            *self._build_options()
        ])

    def _check_command_output(self) -> None:
        """
        Checks if the command executed successfully.
        :return: None
        """
        if not self._command.returncode == 0:
            raise ToolExecutionError(f'Error executing {self.name}: {self.stderr}')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Symlink the input FASTQ files
        FWD_READS = self._tool_inputs['FASTQ_FWD'][0].path
        REV_READS = self._tool_inputs['FASTQ_REV'][0].path

        # Run the command
        self.__build_command(FWD_READS, REV_READS, 'shigatyper_out')
        self._execute_command()

        # Collect the output
        dir_out = self.folder / 'serotype'
        dir_out.mkdir()

        try:
            # Main output
            self._tool_outputs['TSV'] = [ToolIOFile((dir_out / 'shigatyper_out.tsv'))]
            shutil.copy(f'{self.folder}/shigatyper_out.tsv', dir_out)
            # List of hits
            self._tool_outputs['TSV_HITS'] = [ToolIOFile((dir_out / 'shigatyper_out-hits.tsv'))]
            shutil.copy(f'{self.folder}/shigatyper_out-hits.tsv', dir_out)
        except StopIteration:
            raise ToolExecutionError(f"TSV file not found in output folder: {dir_out}")
        self._parse_tsv(self._tool_outputs['TSV'][0].path)

    def _parse_tsv(self, path_tsv: Path) -> None:
        """
        Parses the output TSV file and stores the results in the informs.
        :param path_tsv: Path to output file
        :return: None
        """
        data_serotype = pd.read_table(path_tsv)
        output_dict = data_serotype.to_dict('records')[0]
        self._informs['species'] = output_dict['prediction']
