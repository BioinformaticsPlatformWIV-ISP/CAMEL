from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class ConFindr(Tool):
    """
    This program is designed to find bacterial intra-species contamination in raw Illumina data. It does this by looking
    for multiple alleles of core, single copy genes.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('ConFindr', '0.7.4', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not any(x in self._tool_inputs for x in ('FASTQ_SE', 'FASTQ_PE', 'FASTQ_GZ_SE', 'FASTQ_GZ_PE')):
            raise InvalidInputSpecificationError('FASTQ[_GZ]_(SE/PE) input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        dir_input = self.__symlink_input_files()
        dir_out = Path(self._folder) / 'out'
        self._command.command = ' '.join([
            self._tool_command,
            '--input_directory', str(dir_input),
            '--output_name', str(dir_out),
            '--cross_details',
            '--threads', '1'
        ] + self._build_options())
        self._execute_command()
        self.__set_output(dir_out)

    def __symlink_input_files(self) -> Path:
        """
        Symlinks the input files.
        :return: Path to input directory
        """
        dir_working = Path(self._folder)
        if 'FASTQ_SE' in self._tool_inputs:
            (dir_working / 'sample_R1.fastq').symlink_to(self._tool_inputs['FASTQ_SE'][0].path)
        elif 'FASTQ_GZ_SE' in self._tool_inputs:
            (dir_working / 'sample_R1.fastq.gz').symlink_to(self._tool_inputs['FASTQ_GZ_SE'][0].path)
        elif 'FASTQ_PE' in self._tool_inputs:
            (dir_working / 'sample_R1.fastq').symlink_to(self._tool_inputs['FASTQ_PE'][0].path)
            (dir_working / 'sample_R2.fastq').symlink_to(self._tool_inputs['FASTQ_PE'][1].path)
        else:
            (dir_working / 'sample_R1.fastq.gz').symlink_to(self._tool_inputs['FASTQ_GZ_PE'][0].path)
            (dir_working / 'sample_R2.fastq.gz').symlink_to(self._tool_inputs['FASTQ_GZ_PE'][1].path)
        return dir_working

    def __set_output(self, dir_out: Path) -> None:
        """
        Sets the tool output for this tool.
        :param dir_out: Output directory
        :return: None
        """
        csv_file = dir_out / 'confindr_report.csv'
        self._tool_outputs['CSV'] = [ToolIOFile(csv_file)]
        self._informs.update(pd.read_csv(csv_file).to_dict('records')[0])

    def _check_command_output(self) -> None:
        """
        Checks if the tool was executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Error executing {self.name}: {self._command.stderr}")
