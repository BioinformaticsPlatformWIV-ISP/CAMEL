import os

from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtools import Samtools


class SamtoolsMPileup(Samtools):
    """
    Multi-way pileup.
    Notes:
    - VCF outputs are always bgzipped.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super(SamtoolsMPileup, self).__init__('samtools mpileup', '1.3.1', camel)

    def _check_parameters(self):
        """
        Checks the parameters.
        :return: None
        """
        if self._parameters['output_format'].value not in ['pileup', 'vcf', 'bcf']:
            raise InvalidParameterError("Invalid output format: {}".format(self._parameters['output_format'].value))
        super(SamtoolsMPileup, self)._check_parameters()

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if 'BAM' not in self._tool_inputs:
            raise ValueError("No BAM input file found")
        super(Samtools, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __build_command(self):
        """
        Builds the command.
        :return: None
        """
        command_parts = [
            self._tool_command,
            ' '.join(f.path for f in self._tool_inputs['BAM']),
        ]
        command_parts += self._build_options(['output_format'])
        if 'FASTA' in self._tool_inputs:
            command_parts.append('--fasta-ref {}'.format(self._tool_inputs['FASTA'][0].path))
        if 'TXT_RG' in self._tool_inputs:
            command_parts.append('--exlude-RG {}'.format(self._tool_inputs['TXT_RG'][0].path))
        if 'TXT_POS' in self._tool_inputs:
            command_parts.append('--positions {}'.format(self._tool_inputs['TXT_POS'][0].path))
        if self._parameters['output_format'].value == 'vcf':
            command_parts.append('--VCF')
        elif self._parameters['output_format'].value == 'bcf':
            command_parts.append('--BCF')
        self._command.command = ' '.join(command_parts)

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        output_files = {
            'vcf': ('VCF_GZ', os.path.join(self._folder, self._parameters['output_filename'].value)),
            'bcf': ('BCF', os.path.join(self._folder, self._parameters['output_filename'].value)),
            'pileup': ('PILEUP', os.path.join(self._folder, self._parameters['output_filename'].value))
        }
        key, path = output_files.get(self._parameters['output_format'].value)
        self._tool_outputs[key] = [ToolIOFile(path)]

    def _check_command_output(self):
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if self._command.returncode == 0:
            return
