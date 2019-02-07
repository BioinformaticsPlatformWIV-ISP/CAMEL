import logging
import os

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.components.files.fastautils import FastaUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class ContigQCFilter(Tool):

    """
    Tools that filter contig based on specific criteria.
    - minimal length filter
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: CAMEL instance
        """
        super(ContigQCFilter, self).__init__('ContigQCFilter', '0.1', camel)
        self._contigs = []
        self._contigs_remain = []

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if 'FASTA_contig' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Required FASTA contig input is missing")
        if len(self._tool_inputs['FASTA_contig']) > 1:
            raise InvalidInputSpecificationError("Multiple FASTA contig inputs are invalid, only allow one FASTA file.")

        super(ContigQCFilter, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__readin_contigs()
        self.__filter_contig_by_length()
        self.__set_output()

    def __readin_contigs(self):
        """
        Read in contig sequences as a list
        :return: None
        """
        fasta_file = self._tool_inputs['FASTA_contig'][0].path
        self._contigs = FastaUtils.read_as_dict(fasta_file).values()
        logging.info("{} contigs are read in from FASTA file {!r}.".format(len(self._contigs), os.path.basename(fasta_file)))

    def __filter_contig_by_length(self):
        """
        Filter contigs shorter then 'min_length' parameter
        :return: None
        """
        min_length = int(self._parameters['min_length'].value)
        if min_length == 0:
            logging.info("Filtering contigs by minimal length disabled (min_length cutoff set as 0) ...")
            return
        else:
            logging.info("Filtering contigs by minimal length cutoff ...")
        for contig in self._contigs:
            if len(contig.seq) < min_length:
                logging.info("\tContig {!r} ({}bp) shorter then minimal length required ({}) is removed.".format(contig.id, len(contig.seq), min_length))
            else:
                self._contigs_remain.append(contig)
        logging.info("{} contigs remain after filtering by minimal length cutoff.".format(len(self._contigs_remain)))

        self._contigs = self._contigs_remain
        self._contigs_remain = []

    def __set_output(self):
        """
        Output self._contigs_remain to output FASTA file and set output specs:
        :retrun: None
        """
        fasta_file_out = os.path.join(self._folder, self._parameters['fasta_output'].value)
        FastaUtils.write(self._contigs, fasta_file_out)
        self._tool_outputs['FASTA_contig'] = [ToolIOFile(fasta_file_out)]
