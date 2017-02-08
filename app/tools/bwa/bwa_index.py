import os

from app.components.files.fileutils import FileUtils
from app.io.tooliovalue import ToolIOValue
from app.tools.bwa.bwa import BWA


class BWAIndex(BWA):
    """BWAIndex genome using 'bwa index' from BWA for read mapping"""

    MULTI_FASTA_GENOME_FILE = 'complete_genome.fasta'

    def __init__(self, camel):
        """
        Initialize BWA index
        :param camel: Camel instance
        :return: None
        """
        super(BWAIndex, self).__init__('bwa_index', '0.7.15', camel)
        self._refgenome_fasta = None

    def _execute_tool(self):
        """
        Function to run BWA index
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __get_multi_fasta_genome_filename(self):
        """
        Get the filename used for multi fasta file representing complete genome
        :return: name of the multi fasta file with complete path
        """
        return os.path.join(self._folder, BWAIndex.MULTI_FASTA_GENOME_FILE)

    def _check_input(self):
        """
        Check FASTA_REF input and concatenate them if multiple fasta input files
        :return: None
        """
        super(BWAIndex, self)._check_input()

        nb_of_inputs = len(self._tool_inputs['FASTA_REF'])
        if nb_of_inputs == 0:
            raise ValueError("Required reference genome (FASTA) input file is missing.")
        elif nb_of_inputs > 1:
            multifasta_file = self.__get_multi_fasta_genome_filename()
            FileUtils.concatenate_files(multifasta_file, [f.path for f in self._tool_inputs['FASTA_REF']])
            self._refgenome_fasta = multifasta_file
        else:
            fasta_file_path = self._tool_inputs['FASTA_REF'][0].path
            self._refgenome_fasta = os.path.join(self._folder, os.path.basename(fasta_file_path))
            os.symlink(fasta_file_path, self._refgenome_fasta)

    def __set_output(self):
        """
        Set output for BWA index
        :return: None
        """
        self._tool_outputs['INDEX_GENOME_PREFIX'] = [ToolIOValue(self._refgenome_fasta)]

    def __build_command(self):
        """
        Build the command to run BWA index
        :return: None
        """
        self._command.command = '{} {} {}'.format(
            self._tool_command, ' '.join(self._build_options()), self._refgenome_fasta
        )
