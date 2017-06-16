import os
import re

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.invalidparametererror import InvalidParameterError
from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue
from app.tools.bedtools.bedtoolsgetfasta import BedtoolsGetFasta
from app.components.blasthit.blastnfmt6tsvparser import BlastnFmt6TSVParser


class BedtoolsGetBlastFasta(BedtoolsGetFasta):

    """
    Use Bedtools getfasta function to extract sequences based on BLAST aglinment, targeted sequences can be extracted
    from either BLAST query or subject.
    """
    BED_FILE = 'blasthits_sequences.bed'

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(BedtoolsGetBlastFasta, self).__init__(camel, 'bedtools getblastfasta', '2.25.0')

        self._required_inputs = ['TSV_BLAST', 'FASTA']

        self._specific_parameters = ['mode', 'target', 'outputfile_name']
        self._BED_files = []
        self._targets = []
        self._target_file_prefix = []
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
        for bed_file, fasta_file in zip(self._BED_files, self._FASTA_files):
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
        Check input for bedtools getfasta
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
        Set proper input for bedtools getfasta
        :return: None
        """
        self.__output_blasthits_to_bed(self._tool_inputs['TSV_BLAST'][0].path)

        self._input_specs = ["-fi {}".format(self._tool_inputs['FASTA'][0].path)]

    @staticmethod
    def __retrieve_query_sequence(hit):
        """
        Extract blast query sequences information for extraction
        :param hit: one blast hit information parsed from blast tabular output
        :return: bed information properly formatted
        """
        # customized blast outfmt 6 data columns
        # 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore strand qcovs'
        # BED format: chrom, chromStart, chromEnd, name, score, strand
        if hit.sstart < hit.send:
            bed_info = [hit.qseqid, hit.qstart - 1, hit.qend, hit.sseqid, 0, '+']
        else:
            bed_info = [hit.qseqid, hit.qstart - 1, hit.qend, hit.sseqid, 0, '-']
        return bed_info

    @staticmethod
    def __retrieve_subject_sequence(hit):
        """
        Extract blast subject sequences information for extraction
        :param hit: one blast hit information parsed from blast tablular output
        :return: bed information properly formatted
        """
        # customized blast outfmt 6 data columns
        # 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore strand qcovs'
        # BED format: chrom, chromStart, chromEnd, name, score, strand
        if hit.sstart < hit.send:
            bed_info = [hit.sseqid, hit.sstart - 1, hit.send, hit.qseqid, 0, '+']
        else:
            bed_info = [hit.sseqid, hit.send - 1, hit.sstart, hit.qseqid, 0, '-']
        return bed_info

    @staticmethod
    def __remove_special_characters(string_in):
        """
        Remove special characters in a string so it can be used in filename
        :param string_in: input string might contain special characters
        :return: string with special characters replaced
        """
        string_out = re.sub('[^A-Za-z0-9]+', '_', string_in).strip()
        return string_out

    def __output_blasthits_to_separate_bed(self, blasthits):
        """
        Output blast hits into separate bed files
        :param blasthits: list of blasthits extracted from blast output
        :return: None
        """
        opened_files = {}
        target = self._parameters['target'].value

        try:
            for hit in blasthits:
                if target == 'subject':
                    key = self.__remove_special_characters(hit.qseqid)  # 'query id' as key
                    bed_info = self.__retrieve_subject_sequence(hit)
                elif target == 'query':
                    key = self.__remove_special_characters(hit.sseqid)  # 'subject id'as key
                    bed_info = self.__retrieve_query_sequence(hit)
                else:
                    raise InvalidParameterError(
                        "Unsupported target {!r} for Bedtools GetBlastFasta, should be 'subject' or 'query'.".format(target))

                if key not in opened_files:
                    bed_file = os.path.join(self._folder, key + '-' + self.BED_FILE)
                    self._targets.append(hit.qseqid)
                    self._target_file_prefix.append(key)
                    self._BED_files.append(bed_file)
                    opened_files[key] = open(bed_file, 'w')

                opened_files[key].write("\t".join(map(str, bed_info)) + "\n")

        finally:
            for k, f in opened_files.iteritems():
                f.close()

    def __output_blasthits_to_one_bed(self, blasthits):
        """
        Output blast hits into one bed file
        :param blasthits: list of blasthits extracted from blast output
        :return: None
        """
        bed_file = os.path.join(self._folder, self.BED_FILE)
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
        self._target_file_prefix = []
        self._BED_files.append(bed_file)

    def __output_blasthits_to_bed(self, blasthits_file):
        """
        Output hits dict into one or multiple bed files based on 'mode'
        :param blasthits_file: the tsv file of blast outfmt 6
        :return: None
        """
        blastn_file = BlastnFmt6TSVParser(
            blasthits_file, with_seq=False, columns="qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore sstrand qcovs qcovhsp".split(" "))
        blasthits = blastn_file.read_hits_as_list()
        mode = self._parameters['mode'].value
        if mode == 'all':
            self.__output_blasthits_to_one_bed(blasthits)
        elif mode == 'individual':
            self.__output_blasthits_to_separate_bed(blasthits)
        else:
            raise InvalidParameterError(
                "Unsupported mode {!r} for Bedtools GetBlastFasta, should be 'all' or 'individual'.".format(mode))

    def __set_output(self):
        """
        Set the output specification
        :return: None
        """
        self._FASTA_files = [os.path.splitext(f)[0] + ".fa" for f in self._BED_files]
        self._tool_outputs.update({
            'BED': [ToolIOFile(x) for x in self._BED_files],
            'FASTA': [ToolIOFile(x) for x in self._FASTA_files],
            'Targets': [ToolIOValue(x) for x in self._targets]
        })
