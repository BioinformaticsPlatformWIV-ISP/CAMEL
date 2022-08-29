from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class ResFinder(Tool):
    """
    ResFinder identifies acquired antimicrobial resistance genes in total or partial sequenced isolates of bacteria.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool
        :param camel: CAMEL instance
        """
        super().__init__('ResFinder', '4.1.11', camel)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if any(key in self._tool_inputs for key in ('FASTA' or 'FASTQ_PE' or 'FASTQ_SE')):
            raise InvalidInputSpecificationError('FASTA or FASTQ input is required')
        if not ('acquired' or 'point' in self._tool_inputs):
            raise InvalidInputSpecificationError('Either "acquired" or "point" is required')
        super()._check_input()

    def _build_command(self) -> None:
        """
        Builds the command to run resfinder.
        :return: None
        """
        if 'FASTA' in self._tool_inputs:
            input_str = f'--inputfasta {self._tool_inputs["FASTA"][0].path}'
        elif 'FASTQ_SE' in self._tool_inputs:
            input_str = f'--inputfastq {self._tool_inputs["FASTQ"][0].path}'
        elif 'FASTQ_PE' in self._tool_inputs:
            input_str = f'--inputfastq {self._tool_inputs["FASTQ_PE"][0].path} ' \
             f'{str(self._tool_inputs["FASTQ_PE"][1].path)}'
        self._command.command = ' '.join([self._tool_command, input_str, ' '.join(self._build_options())])

    def _check_command_output(self) -> None:
        """
        Checks command output.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    def _set_output(self) -> None:
        """
        Collects the tool output.
        """
        dir_out = self.folder / self._parameters['output_path'].value

        self._tool_outputs['TSV_pheno_general'] = [
            ToolIOFile(dir_out / Path('pheno_table.txt'))
        ]
        if 'acquired' in self._parameters:
            self._tool_outputs['TSV_genes'] = [
                ToolIOFile(dir_out / Path('ResFinder_results_tab.txt'))
            ]
        if 'point' in self._parameters:
            self._tool_outputs['TSV_point'] = [
                ToolIOFile(dir_out / Path('PointFinder_results.txt'))
            ]
            specific_species = '_'.join(self._parameters['species'].value.lower().split(' ')).replace('"', '')
            self._tool_outputs['TSV_pheno_species'] = [
                ToolIOFile(dir_out / Path(f'pheno_table_{specific_species}.txt'))
            ]

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._set_output()
