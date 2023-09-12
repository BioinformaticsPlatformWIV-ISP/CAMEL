from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Clair3(Tool):
    """
    Clair3 is a germline small variant caller for long-reads.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool
        :param camel: CAMEL instance
        """
        super().__init__('Clair3', '1.0.4', camel)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTA reference is required')
        if 'BAM' not in self._tool_inputs:
            raise InvalidInputSpecificationError('BAM alignment file is required')

        input_folder = self._tool_inputs['FASTA'][0].path.parent
        base_fasta_name = self._tool_inputs['FASTA'][0].path.name
        fasta_index_file = [f for f in input_folder.glob(f'{base_fasta_name}.fai')]
        if not (len(fasta_index_file) > 0):
            raise InvalidInputSpecificationError('FASTA reference needs to be indexed')
        super()._check_input()

    def _build_command(self, fasta_input: Path, bam_input: Path) -> None:
        """
        Builds the command to run clair3.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command, f'--bam_fn {bam_input}', f'--ref_fn {fasta_input}',
            *self._build_options(delimiter='=',
                                 excluded_parameters=['haploid_precise', 'no_phasing', 'include_ctgs', 'long_indel']),
            *self._build_options(excluded_parameters=['chunk_size', 'model_path', 'platform', 'threads', 'output_path'])
        ])

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
        self._tool_outputs['VCF'] = [ToolIOFile(dir_out / Path('merge_output.vcf.gz'))]

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        fasta_input = Path(str(self._tool_inputs['FASTA'][0]))
        bam_input = Path(str(self._tool_inputs['BAM'][0]))
        self._build_command(fasta_input, bam_input)
        self._execute_command()
        self._set_output()
