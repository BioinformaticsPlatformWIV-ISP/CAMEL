import os

from app.components.files.fileutils import FileUtils
from app.io.tooliovalue import ToolIOValue
from app.tools.bowtie2.bowtie2 import Bowtie2


class Index(Bowtie2):
    """
    Index genome using 'bowtie2-build' cmd of Bowtie2
    """

    MULTI_FASTA_GENOME_FILE = 'complete_genome.fasta'

    def __init__(self, camel):
        """
        Initialize bowtie2 index
        :param camel: Camel instance
        :return: none
        """
        super(Index, self).__init__('bowtie2 index', '2.3.0', camel)
        self._refgenome_fasta = None

    def _execute_tool(self):
        """
        Function to run BWA index
        :return: none
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __get_multi_fasta_genome_filename(self):
        """
        Get the filename used for multi fasta file representing complete genome
        :return: name of the multi fasta file with complete path
        """
        return os.path.join(self._folder, self.MULTI_FASTA_GENOME_FILE)

    def _check_input(self):
        """
        Check FASTA_REF input and concatenate them if multiple fasta input files
        :return: none
        """
        no_of_inputs = len(self._tool_inputs['FASTA_REF'])
        if no_of_inputs > 1:
            multifasta_file = self.__get_multi_fasta_genome_filename()
            FileUtils.concatenate_files(multifasta_file, [f.path for f in self._tool_inputs['FASTA_REF']])
            self._refgenome_fasta = multifasta_file
        elif no_of_inputs == 0:
            raise ValueError("Required reference genome (FASTA) input file is missing.")
        else:
            self._refgenome_fasta = self._tool_inputs['FASTA_REF'][0].path

    def __set_output(self):
        """
        Set output for bowtie2 index
        :return: none
        """
        self._tool_outputs['INDEX_GENOME_PREFIX'] = [ToolIOValue(self._refgenome_fasta)]

    def __build_command(self):
        """
        Build the command to run bowtie2 index
        :return: none
        """
        # Note the refgenome fasta name is used as index base
        self._command.command = "{} {} {} {}".format(
            self._tool_command,
            " ".join(self._build_options()),
            self._refgenome_fasta,
            self._refgenome_fasta
        )
