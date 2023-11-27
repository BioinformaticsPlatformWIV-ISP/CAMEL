import os
import re

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.bedtools.bedtoolsgetfasta import BedtoolsGetFasta


class BedtoolsGetBlastFasta(BedtoolsGetFasta):

    """
    Use Bedtools getfasta function to extract sequences based on BLAST alignment (outfmt6 TSV) and query/subject
    sequences (FASTA). The tool can extract all hit sequences into one fasta file, or alternatively, group hit
    sequences per target (BLAST query or subject) and output into separate fasta files. Note that query sequence MUST
    be extracted from query fasta file and subject sequence should be extracted from subject fasta file.
    """
    # TODO: update blastnhit related code upon update blastnfmt6tsvparser
    DEFAULT_BEDFILE_NAME = 'blasthits_sequences.bed'

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__(camel, 'bedtools getblastfasta', '2.25.0')

        self._required_inputs = ['TSV_BLAST', 'FASTA']

        self._specific_parameters = ['mode', 'target', 'outputfile_name']
        self._bed_files = []
        self._targets = []
        self._target_file_prefixes = []
        self._input_specs = []
        self._FASTA_files = []
        self._input_string = ''
        self._output_string = ''

    def _execute_tool(self):
        """
        Executes this tool
        :return: None
        """
        self.__set_input()
        self.__set_output()
        for bed_file, fasta_file in zip(self._bed_files, self._FASTA_files):
            self.__run_getfasta(bed_file, fasta_file)

    def __run_getfasta(self, bed_file, fasta_file):
        """
        Run Bedtools getfasta command to extract sequences based on input BED files
        :return: None
        """
        self._input_string = " ".join(self._input_specs) + " -bed {}".format(bed_file)
        self._output_string = "-fo {}".format(fasta_file)
        self.__build_command()
        self._execute_command()

    def __build_command(self):
        """
        Build the command to run tool
        :return: None
        """
        self._command.command = " ".join([
            self._tool_command,
            self._input_string,
            self._output_string,
            " ".join(self._build_options(excluded_parameters=self._specific_parameters))
        ])

    def _check_input(self):
        """
        Check the input specifications
        :return: None
        """
        self._check_required_inputs()

        if len(self._tool_inputs['TSV_BLAST']) != 1:
            raise InvalidInputSpecificationError("Exactly one TSV_BLAST input file expected.")
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidInputSpecificationError("Exactly one FASTA input file expected.")

        super(BedtoolsGetFasta, self)._check_input()

    def __set_input(self):
        """
        Set input(s) for bedtools getfasta
        :return: None
        """
        self.__output_blasthits_to_bed(self._tool_inputs['TSV_BLAST'][0].path)

        self._input_specs = ["-fi {}".format(self._tool_inputs['FASTA'][0].path)]

    def __set_output(self):
        """
        Set the output specifications
        :return: None
        """
        self._FASTA_files = [os.path.splitext(f)[0] + ".fa" for f in self._bed_files]
        self._tool_outputs.update({
            'BED': [ToolIOFile(x) for x in self._bed_files],
            'FASTA': [ToolIOFile(x) for x in self._FASTA_files],
            'Targets': [ToolIOValue(x) for x in self._targets]
        })

    def __output_blasthits_to_bed(self, blasthits_file):
        """
        Output hits dict into one or multiple bed files based on 'mode'
        :param blasthits_file: blastn outfmt 6 tsv file
        :return: None
        """
        blastn_file = BlastnFmt6TSVParser(
            blasthits_file, columns="qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore sstrand qcovs qcovhsp".split(" "))
        blasthits = blastn_file.read_hits_as_list()
        mode = self._parameters['mode'].value
        if mode == 'all':
            self.__output_blasthits_to_one_bed(blasthits)
        elif mode == 'individual':
            self.__output_blasthits_to_separate_bed(blasthits)
        else:
            raise InvalidParameterError(
                "Unsupported mode {!r} for Bedtools GetBlastFasta, should be 'all' or 'individual'.".format(mode))

    def __output_blasthits_to_separate_bed(self, blasthits):
        """
        Output blast hits into separate bed files (per target)
        :param blasthits: list of blasthits extracted from blast output
        :return: None
        """
        opened_files = {}
        target = self._parameters['target'].value

        try:
            for hit in blasthits:
                if target == 'subject':
                    target_id = hit.qseqid
                    key = self.__remove_special_characters(hit.qseqid)  # 'query id' as key
                    bed_info = self.__retrieve_subject_sequence(hit)
                elif target == 'query':
                    target_id = hit.qseqid
                    key = self.__remove_special_characters(hit.sseqid)  # 'subject id'as key
                    bed_info = self.__retrieve_query_sequence(hit)
                else:
                    raise InvalidParameterError(
                        "Unsupported target {!r} for Bedtools GetBlastFasta, should be 'subject' or 'query'.".format(target))

                if key not in opened_files:
                    bed_file = os.path.join(self._folder, key + '-' + self.DEFAULT_BEDFILE_NAME)
                    self._targets.append(target_id)
                    self._target_file_prefixes.append(key)
                    self._bed_files.append(bed_file)
                    opened_files[key] = open(bed_file, 'w')

                opened_files[key].write("\t".join(map(str, bed_info)) + "\n")

        finally:
            for k, f in opened_files.items():
                f.close()

    def __output_blasthits_to_one_bed(self, blasthits):
        """
        Output blast hits into one bed file
        :param blasthits: list of blasthits extracted from blast output
        :return: None
        """
        bed_file = os.path.join(self._folder, self.DEFAULT_BEDFILE_NAME)
        target = self._parameters['target'].value

        with open(bed_file, 'w') as outf:
            for hit in blasthits:
                if target == 'subject':
                    bed_info = self.__retrieve_subject_sequence(hit)
                elif target == 'query':
                    bed_info = self.__retrieve_query_sequence(hit)
                else:
                    raise InvalidParameterError(
                        "Unsupported target {!r} for Bedtools GetBlastFasta, should be 'subject' or 'query'.".format(target))

                outf.write("\t".join(map(str, bed_info)) + "\n")

        self._targets = []
        self._target_file_prefixes = []
        self._bed_files.append(bed_file)

    @staticmethod
    def __remove_special_characters(string_in):
        """
        Remove special characters in a string so it can be used in filename
        :param string_in: input string
        :return: string with special characters replaced
        """
        return re.sub('[^A-Za-z0-9]+', '_', string_in).strip()

    @staticmethod
    def __retrieve_query_sequence(hit):
        """
        Extract blastn query hit sequences
        :param hit: one blast hit information parsed from blast tabular output
        :return: bed information properly formatted
        """
        # blastn hits with qseqid, sseqid, qstart, qend, sstart, send
        #
        # BED format: chrom, chromStart, chromEnd, name, score, strand
        return [hit.qseqid, hit.qstart - 1, hit.qend, hit.sseqid, 0, '+']

    @staticmethod
    def __retrieve_subject_sequence(hit):
        """
        Extract blastn subject hit sequences
        :param hit: one blast hit information parsed from blast tablular output
        :return: bed information properly formatted
        """
        # blastn hits with qseqid, sseqid, qstart, qend, sstart, send
        #
        # BED format: chrom, chromStart, chromEnd, name, score, strand
        if hit.sstart < hit.send:
            bed_info = [hit.sseqid, hit.sstart - 1, hit.send, hit.qseqid, 0, '+']
        else:
            bed_info = [hit.sseqid, hit.send - 1, hit.sstart, hit.qseqid, 0, '-']
        return bed_info
