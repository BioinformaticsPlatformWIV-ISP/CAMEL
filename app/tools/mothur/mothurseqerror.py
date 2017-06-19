from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliofile import ToolIOFile
from app.tools.mothur.mothur import Mothur


class MothurSeqError(Mothur):
    """
    Seq.error measures the error rate
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(MothurSeqError, self).__init__('mothur_seq_error', '1.39.1', camel)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA and FASTA_REF are required
        - For now 'TSV_Names', 'TSV_Counts', and 'TSV_Qfile' are not implemented (lack of documentation)
        - Only one file per key allowed
        :return: None
        """
        super(MothurSeqError, self)._check_input()
        if 'FASTA' not in self._tool_inputs or 'FASTA_REF' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Invalid input files (keys) given for Mothur '
                                                 'seq.error: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.iteritems():
            if key in ['TSV_Names', 'TSV_Counts', 'TSV_Qfile']:
                raise InvalidInputSpecificationError('These input keys are not yet implemented as documentation for these keys'
                                                     ' is not yet available: {!r}'.format(self._tool_inputs))
            elif key not in ['FASTA', 'FASTA_REF']:
                raise InvalidInputSpecificationError('Invalid input key given for Mothur seq.error: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidInputSpecificationError('Invalid number (max = 1) of files in each key given for Mothur \
                                                     seq.error: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = ['fasta={}'.format(self._tool_inputs['FASTA'][0]),
                 'reference={}'.format(self._tool_inputs['FASTA_REF'][0]),
                 'outputdir={}'.format(self._folder)]
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = super(MothurSeqError, self)._get_basename()
        self._tool_outputs.update({
            'TSV_Summary': [ToolIOFile(basename + '.error.summary')],
            'TSV_Seq': [ToolIOFile(basename + '.error.seq')],
            'TSV_Chimera': [ToolIOFile(basename + '.error.chimera')],
            'TSV_SeqForward': [ToolIOFile(basename + '.error.seq.forward')],
            'TSV_SeqReverse': [ToolIOFile(basename + '.error.seq.reverse')],
            'TSV_Count': [ToolIOFile(basename + '.error.count')],
            'TSV_Matrix': [ToolIOFile(basename + '.error.matrix')]
        })
