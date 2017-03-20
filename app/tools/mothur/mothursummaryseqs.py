from app.tools.mothur.mothur import Mothur
from app.io.tooliofile import ToolIOFile
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError


class MothurSummarySeqs(Mothur):
    """
    The summary.seqs command will summarize the quality of sequences in
    an unaligned or aligned fasta-formatted sequence file.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(MothurSummarySeqs, self).__init__('mothur_summary_seqs', '1.39.1', camel)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - Allowed keys are 'FASTA', 'TSV_Counts'
        - Only one input file allowed per key
        :return: None
        """
        super(MothurSummarySeqs, self)._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('No input file given for Mothur summary.seqs: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.iteritems():
            if key not in ['FASTA', 'TSV_Counts']:
                raise InvalidInputSpecificationError('Invalid input key given for Mothur summary.seqs: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidInputSpecificationError('Invalid number (max = 1) of files given for Mothur \
                                                     summary.seqs: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and input/output directories
        Example: fasta=File1.trim.contig.fasta, inputdir=/test/data/input/,
        outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        items = ['fasta={}'.format(self._tool_inputs['FASTA'][0])]
        if 'TSV_Counts' in self._tool_inputs:
            items.append('count={}'.format(self._tool_inputs['TSV_Counts'][0]))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the output file object with them
        :return: None
        """
        basename = super(MothurSummarySeqs, self)._get_basename()
        self._tool_outputs['TSV_Summary'] = [ToolIOFile(basename + '.summary')]
        self._tool_outputs['TSV_Stats'] = [ToolIOFile(basename + '.stats')]

    def _execute_tool(self):
        """
        Runs Prinseq
        :return: None
        """
        self._create_symlinks()
        self._build_command()
        self._execute_command()
        self.__write_stats_to_file()
        self._remove_symlinks()
        self._set_output()

    def __write_stats_to_file(self):
        """
        Writes the statistics that were output to stdout to a file
        :return: None
        """
        output_file = file(self._get_basename() + '.stats', 'w')
        write_to_output = False
        for line in self._command.stdout.splitlines():
            # The first line of the stats starts with two tabs (i.e. \t\tStart\tEnd...)
            if line.startswith('\t'):
                write_to_output = True
            if write_to_output is True:
                output_file.write(line + '\n')
            # The last line of the stats starts with a '#' (i.e. # of Seqs...)
            if line.startswith('#'):
                break
