from app.tools.mothur.mothur import Mothur
from app.io.tooliofile import ToolIOFile
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError


class MothurClassifyOtu(Mothur):
    """
    The classify.otu command is used to get a consensus taxonomy for an otu.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(MothurClassifyOtu, self).__init__('mothur_classify_otu', '1.39.1', camel)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - TSV_List and TSV_Taxonomy keys are required
        - Additional allowed keys: 'TSV_Groups', 'TSV_Counts', 'TSV_Names', 'TSV_RefTaxonomy'
        - Only one input file per key allowed
        :return: None
        """
        super(MothurClassifyOtu, self)._check_input()
        if 'TSV_List' not in self._tool_inputs or 'TSV_Taxonomy' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Invalid input files (keys) given for Mothur '
                                                 'classify.otu: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.iteritems():
            if key not in ['TSV_List', 'TSV_Taxonomy', 'TSV_Groups', 'TSV_Counts', 'TSV_Names', 'TSV_RefTaxonomy']:
                raise InvalidInputSpecificationError('Invalid input key given for Mothur '
                                                     'classify.otu: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidInputSpecificationError('Invalid number (max = 1) of files given for Mothur \
                                                     classify.otu: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = ['list={}'.format(self._tool_inputs['TSV_List'][0]),
                 'taxonomy={}'.format(self._tool_inputs['TSV_Taxonomy'][0])]
        if 'TSV_Counts' in self._tool_inputs:
            items.append('count={}'.format(self._tool_inputs['TSV_Counts'][0]))
        if 'TSV_Groups' in self._tool_inputs:
            items.append('group={}'.format(self._tool_inputs['TSV_Groups'][0]))
        if 'TSV_Names' in self._tool_inputs:
            items.append('name={}'.format(self._tool_inputs['TSV_Names'][0]))
        if 'TSV_RefTaxonomy' in self._tool_inputs:
            items.append('reftaxonomy={}'.format(self._tool_inputs['TSV_RefTaxonomy'][0]))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them.
        :return: None
        """
        labels = super(MothurClassifyOtu, self)._get_labels()
        self._tool_outputs.update({'TSV_Taxonomy': [], 'TSV_Summary': []})
        basename = super(MothurClassifyOtu, self)._get_basename('TSV_List')
        # Each label creates a seperate output
        for label in labels:
            self._tool_outputs['TSV_Taxonomy'] += [ToolIOFile(basename + '.' + label + '.cons.taxonomy')]
            self._tool_outputs['TSV_Summary'] += [ToolIOFile(basename + '.' + label + '.cons.tax.summary')]
