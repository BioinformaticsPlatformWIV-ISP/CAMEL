import logging
import os

from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue
from app.tools.tool import Tool


class Srst2Mlst(Tool):
    """
    This program is designed to take Illumina sequence data, a MLST database and/or a database of gene sequences
    (e.g. resistance genes, virulence genes, etc) and report the presence of subtypes and/or reference genes.
    """

    def __init__(self, camel):
        """
        Initialize SRST2 MLST tool.
        :param camel: Camel instance
        """
        super(Srst2Mlst, self).__init__('SRST2_MLST', '0.1.5', camel)

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self._command.command = self.__build_command()
        self._execute_command()
        self.__set_output()

    def __build_command(self):
        """
        Builds the command line command.
        :return: Command line command
        """
        return ' '.join([self._tool_command,
                         self.__build_input_string(),
                         '--mlst_db {}'.format(self._tool_inputs['FASTA'][0].path),
                         ' '.join(self._build_options(excluded_parameters=['threads'])),
                         '--other "--threads {}"'.format(self._parameters['threads'].value)])

    def __build_input_string(self):
        """
        Builds a string containing the input.
        :return: Input options string
        """
        if 'FASTQ_PE' in self._tool_inputs:
            return '--input_pe {}'.format(' '.join([f.path for f in self._tool_inputs['FASTQ_PE']]))
        else:
            return '--input_se {}'.format(self._tool_inputs['FASTQ_SE'][0].path)

    def __set_output(self):
        """
        Sets the output files.
        :return: None
        """
        for file_ in os.listdir(self._folder):
            full_path = os.path.join(self._folder, file_)
            key = self._get_output_file_key(file_)
            if key is not None:
                self._tool_outputs[key] = [ToolIOFile(full_path)]
        self._tool_outputs['VAL_Sequence_type'] = [
            ToolIOValue(Srst2Mlst.__get_sequence_type(self._tool_outputs['TSV'][0].path))]

    def _check_input(self):
        """
        Checks whether the given inputs are valid.
        - FASTQ_PE or FASTQ_SE reads are required (checked by super class)
        - FASTA file with allele sequences is required
        - MLST file with sequence type definitions is optional
        """
        super(Srst2Mlst, self)._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise IOError('No FASTA file with MLST alleles found.')
        if 'MLST' not in self._tool_inputs:
            logging.info("No MLST definitions found. Only performing allele detection.")

    def _get_output_file_key(self, filename):
        """
        Returns the key for the given output file.
        :param filename: Filename
        :return: Key
        """
        output_filename = self._parameters['output_filename'].value
        if all([x in filename for x in ['mlst', 'results']]):
            return 'TSV'
        elif filename.endswith('.pileup') and filename.startswith(output_filename):
            return 'PILEUP'
        elif filename.endswith('.bam'):
            return 'BAM'
        elif filename.endswith('.scores'):
            return 'TSV_Scores'
        elif filename.endswith('consensus_alleles.fasta'):
            return 'FASTA'
        elif filename.endswith('.log'):
            return 'LOG'

    @staticmethod
    def __get_sequence_type(output_file):
        """
        Parses the output file to obtain the sequence type.
        :param output_file: Output file
        :return: Sequence type
        """
        with open(output_file) as f:
            content = f.readlines()
            if len(content) == 1:
                return 'ND'
            elif len(content) == 2:
                return content[1].split('\t')[1]
            raise ValueError("Invalid SRST2 output file. Content: '{}'".format(content))
