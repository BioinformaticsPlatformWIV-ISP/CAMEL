import logging
import os

from app.components.files.fileutils import FileUtils
from app.io.tooliofile import ToolIOFile
from app.tools.seqtk.seqtk import Seqtk


class SeqtkSubsample(Seqtk):

    """
    Class that subsamples fastq/fasta file(s) using seqkt
    """

    def __init__(self, camel):
        """
        Initialize seqtk
        :param camel: Camel instance
        :return: None
        """
        super(SeqtkSubsample, self).__init__('Seqtk Subsample', '1.2', camel)
        self._supported_inputs = ['FASTA', 'FASTQ', 'FASTA_PE', 'FASTQ_PE']
        self._function_name = 'Subsample'
        self._specific_parameters = ['combine_output', 'output_prefix', 'fraction']
        self._output_files = []

    def _execute_tool(self):
        """
        Function to run seqtk subsample
        :return: None
        """
        self.__set_cmd_output()
        logging.debug("Seqtk Subsample input informs: input_mode {}, input_file_type {}".format(
            self.input_mode, self.input_file_type))
        for idx, value in enumerate(self._input_files):
            self.__build_command_with_iofiles(self._input_files[idx], self._output_files[idx])
            self._execute_command()
        self._set_output()

    def __set_cmd_output(self):
        """
        Set the output specification for run seqtk subsample command
        :return: None
        """
        output_file_prefix = self._parameters['output_prefix'].value
        output_suffix = self.__get_output_file_suffix()

        if self.input_mode == 'SE':
            self._output_files = [os.path.join(self._folder, output_file_prefix + output_suffix)]
        elif self.input_mode == 'PE':
            self._output_files = [
                os.path.join(self._folder, output_file_prefix + "_1" + output_suffix),
                os.path.join(self._folder, output_file_prefix + "_2" + output_suffix)
            ]

    def __build_command_with_iofiles(self, input_file, output_file):
        """
        Build the command to seqtk subsample
        :param input_file:
        :param output_file:
        :return: None
        """
        self._command.command = "{} {} {} {} > {}".format(
            self._tool_command,
            " ".join(self._build_options(excluded_parameters=self._specific_parameters)),
            input_file,
            self._parameters['fraction'].value,
            output_file
        )

    def _set_output(self):
        """
        Set self._tool_outputs specification
        """
        if self.input_mode == 'PE':
            self.__set_pe_outputs()
        else:
            self.__set_se_output()

    def __set_se_output(self):
        """
        Set self._tool_outputs specification for SE mode
        """
        self._tool_outputs[self.input_file_type] = [ToolIOFile(self._output_files[0])]

    def __set_pe_outputs(self):
        """
        Set self._tool_output for PE mode
        :return: None
        """
        # for PE mode, final output depend on combine_output option
        if 'combine_output' in self._parameters:
            # Tag if set: combine outputs of individual PE reads into one file
            output_file = os.path.join(self._folder,
                                       self._parameters['output_prefix'].value + self.__get_output_file_suffix())
            FileUtils.concatenate_files(output_file, self._output_files)
            self._tool_outputs[self.input_file_type] = [ToolIOFile(output_file)]

        else:
            self._tool_outputs[self.input_file_type + '_PE'] = [ToolIOFile(f) for f in self._output_files]

    def __get_output_file_suffix(self):
        """
        Obtain the file's suffix based on input file type
        :return: string, filename suffix
        """
        if self.input_file_type == 'FASTA':
            return '.fa'
        elif self.input_file_type == 'FASTQ':
            return '.fq'
        else:
            raise ValueError("Seqtk Subsample supports only FASTA/FASTQ as input.")
