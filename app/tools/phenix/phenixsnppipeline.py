import os

from app.components.filesystemhelper import FileSystemHelper
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class PhenixSnpPipeline(Tool):
    """
    Public Health England SNP calling pipeline.
    """

    DEFAULT_SAMPLE_NAME = 'test_sample'

    def __init__(self, camel):
        """
        Initializes the SNP calling pipeline.
        :param camel: CAMEL instance
        """
        super(PhenixSnpPipeline, self).__init__('PHEnix SNP Pipeline', '1.2', camel)
        self._sample_name = None

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No reference FASTA input found.")
        if 'FASTQ_PE' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No PE FASTQ input dataset found.")
        if len(self._tool_inputs['FASTQ_PE']) != 2:
            raise InvalidInputSpecificationError("Two FASTQ files (forward & reverse) are expected.")
        super(PhenixSnpPipeline, self)._check_input()

    def __set_sample_name(self):
        """
        Sets the sample name.
        :return: None
        """
        if 'SAMPLE_NAME' not in self._tool_inputs:
            self._sample_name = PhenixSnpPipeline.DEFAULT_SAMPLE_NAME
        else:
            self._sample_name = FileSystemHelper.make_valid(self._tool_inputs['SAMPLE_NAME'][0].value)

    def _execute_tool(self):
        """
        Executes the tool.
        :return: None
        """
        self.__set_sample_name()
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __build_command(self):
        """
        Builds the tool command.
        :return: None
        """
        self._command.command = '. $VIRTUALENV; {}'.format(' '.join(
            [self._tool_command,
             '--reference {}'.format(self._tool_inputs['FASTA'][0].path),
             '-r1 {} -r2 {}'.format(*[f.path for f in self._tool_inputs['FASTQ_PE']]),
             '--sample-name "{}"'.format(self._sample_name),
             '--outdir {}'.format(self._folder),
             ' '.join(self._build_options())]))

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        output_mapping = {'{}.vcf': 'VCF', '{}.filtered.vcf': 'VCF_Filt', '{}.bam': 'BAM'}
        for key in output_mapping:
            full_path = os.path.join(self._folder, key.format(self._sample_name))
            if os.path.isfile(full_path):
                self._tool_outputs[output_mapping.get(key)] = [ToolIOFile(full_path)]

    def _check_command_output(self):
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Error executing PHEnix SNP calling pipeline: {}".format(self.stderr))
