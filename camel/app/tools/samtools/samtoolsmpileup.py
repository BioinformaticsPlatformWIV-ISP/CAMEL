from camel.app.camel import Camel
from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtools import Samtools


class SamtoolsMPileup(Samtools):
    """
    Multi-way pileup.
    Notes:
    - VCF outputs are always bgzipped.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('samtools mpileup', '1.9', camel)

    def _check_parameters(self) -> None:
        """
        Checks the parameters.
        :return: None
        """
        if self._parameters['output_format'].value not in ['pileup', 'vcf', 'bcf']:
            raise InvalidParameterError(f"Invalid output format: {self._parameters['output_format'].value}")
        super(SamtoolsMPileup, self)._check_parameters()

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if 'BAM' not in self._tool_inputs:
            raise ValueError("No BAM input file found")
        super(Samtools, self)._check_input()

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
            ' '.join(str(f.path) for f in self._tool_inputs['BAM']),
        ]
        command_parts += self._build_options(['output_format'])
        if 'FASTA' in self._tool_inputs:
            command_parts.append(f'--fasta-ref {self._tool_inputs["FASTA"][0].path}')
        if 'TXT_RG' in self._tool_inputs:
            command_parts.append(f'--exlude-RG {self._tool_inputs["TXT_RG"][0].path}')
        if 'TXT_POS' in self._tool_inputs:
            command_parts.append(f'--positions {self._tool_inputs["TXT_POS"][0].path}')
        if self._parameters['output_format'].value == 'vcf':
            command_parts.append('--VCF')
        elif self._parameters['output_format'].value == 'bcf':
            command_parts.append('--BCF')
        self._command.command = ' '.join(command_parts)

    def __set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        output_files = {
            'vcf': ('VCF_GZ', self.folder / self._parameters['output_filename'].value),
            'bcf': ('BCF', self.folder / self._parameters['output_filename'].value),
            'pileup': ('PILEUP', self.folder / self._parameters['output_filename'].value)
        }
        key, path = output_files.get(self._parameters['output_format'].value)
        self._tool_outputs[key] = [ToolIOFile(path)]

    def _check_command_output(self) -> None:
        """
        Checks the command output.
        Supersedes function in Tool class because warnings printed to stderr can cause false abort.
        """
        self._check_stderr()

        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")
