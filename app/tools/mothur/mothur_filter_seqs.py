from app.tools.mothur.mothur import Mothur
from app.io.tooliofile import ToolIOFile
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError


class MothurFilterSeqs(Mothur):
    """
    filter.seqs removes columns from alignments based on a criteria defined by the user.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(MothurFilterSeqs, self).__init__('mothur_filter_seqs', '1.39.1', camel)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - Only one FASTA file allowed
        - No other input keys allowed
        :return: None
        """
        super(MothurFilterSeqs, self)._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('No valid input files given for Mothur filter.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidInputSpecificationError('Invalid number (max = 1) of files given for Mothur \
                                                 filter.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidInputSpecificationError('Too many input keys given voor Mothur filter.seqs: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        input_string = 'fasta={}'.format(self._tool_inputs['FASTA'][0])
        input_string += ', outputdir={}'.format(self._folder)
        return input_string

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = super(MothurFilterSeqs, self)._get_basename()
        self._tool_outputs['FASTA'] = [ToolIOFile(basename + '.filter.fasta')]
