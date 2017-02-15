from app.tools.mothur.mothur import Mothur
from app.io.tooliofile import ToolIOFile
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError


class MothurGetOturep(Mothur):
    """
    The get.oturep command generates a fasta-formatted sequence file containing only a representative sequence
    for each OTU.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(MothurGetOturep, self).__init__('mothur_get_oturep', '1.39.1', camel)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA is required
        - DIST_Phy or DIST_Col key is required
        - TSV_List key is required
        - Additional allowed keys: 'TSV_Counts', 'FASTA', 'TSV_Names', 'TSV_Groups'
        - Only one input file per key allowed
        :return: None
        """
        super(MothurGetOturep, self)._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('No fasta file given for Mothur get.oturep: {!r}'.format(self._tool_inputs))
        if not (('PHY' in self._tool_inputs) != ('DIST' in self._tool_inputs)) \
                or 'TSV_List' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Missing input files (keys) for Mothur '
                                                 'get.oturep: {!r}'.format(self._tool_inputs))
        if 'TSV_Counts' in self._tool_inputs and 'TSV_Names' in self._tool_inputs:
            raise InvalidInputSpecificationError('Both Count and Names file is not allowed for Mothur '
                                                 'get.oturep: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.iteritems():
            if key not in ['PHY', 'DIST', 'FASTA', 'TSV_Counts', 'TSV_Names', 'TSV_Groups', 'TSV_List']:
                raise InvalidInputSpecificationError('Invalid input key given for Mothur get.oturep: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidInputSpecificationError('Invalid number (max = 1) of files given for Mothur \
                                                     get.oturep: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        input_string = ''
        if 'PHY' in self._tool_inputs:
            input_string += 'phylip={}'.format(self._tool_inputs['PHY'][0])
        else:
            input_string += 'column={}'.format(self._tool_inputs['DIST'][0])
        input_string += ', list={}'.format(self._tool_inputs['TSV_List'][0])
        if 'TSV_Names' in self._tool_inputs:
            input_string += ', name={}'.format(self._tool_inputs['TSV_Names'][0])
        elif 'TSV_Counts' in self._tool_inputs:
            input_string += ', count={}'.format(self._tool_inputs['TSV_Counts'][0])
        if 'FASTA' in self._tool_inputs:
            input_string += ', fasta={}'.format(self._tool_inputs['FASTA'][0])
        if 'TSV_Groups' in self._tool_inputs:
            input_string += ', group={}'.format(self._tool_inputs['TSV_Groups'][0])
        input_string += ', outputdir={}'.format(self._folder)
        return input_string

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = super(MothurGetOturep, self)._get_basename('TSV_List')
        self._tool_outputs['FASTA'] = [ToolIOFile(basename + '.unique.rep.fasta')]
        labels = super(MothurGetOturep, self)._get_labels()
        if 'TSV_Counts' in self._tool_inputs:
            self._tool_outputs['TSV_Counts'] = [ToolIOFile(basename + '.unique.rep.count_table')]
        elif 'TSV_Names' in self._tool_inputs:
            self._tool_outputs['TSV_Names'] = []
            for label in labels:
                self._tool_outputs['TSV_Names'] += [ToolIOFile(basename + '.' + label + '.rep.names')]
