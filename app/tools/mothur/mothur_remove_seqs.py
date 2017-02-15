from app.tools.mothur.mothur import Mothur
from app.io.tooliofile import ToolIOFile
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError


class MothurRemoveSeqs(Mothur):
    """
    The remove.seqs command takes a list of sequence names and either a fastq, fasta, name, group, list, count or
    align.report file to generate a new file that does not contain the sequences in the list.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(MothurRemoveSeqs, self).__init__('mothur_remove_seqs', '1.39.1', camel)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - TSV_Accnos is required
        - In addition only two extra keys are permitted
        - Possible additional key comes from: 'FASTA', 'TSV_Names', 'TSV_Counts', 'TSV_Groups',
          'TSV_AlignReport', 'TSV_List', 'TSV_Taxonomy', 'TSV_Qfile', 'FASTQ'
        - Only one file per key is allowed
        :return: None
        """
        super(MothurRemoveSeqs, self)._check_input()
        if 'TSV_Accnos' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Not enough valid input files given for Mothur '
                                                 'remove.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) > 3:
            raise InvalidInputSpecificationError('Too many input keys given for Mothur remove.seqs: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.iteritems():
            if key not in ['TSV_Accnos', 'FASTA', 'TSV_Names', 'TSV_Counts', 'TSV_Groups',
                           'TSV_AlignReport', 'TSV_List', 'TSV_Taxonomy', 'TSV_Qfile', 'FASTQ']:
                raise InvalidInputSpecificationError('Invalid input key given for Mothur remove.seqs: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidInputSpecificationError('Invalid number (max = 1) of files given for Mothur \
                                                     remove.seqs: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        # TSV_Accnos is required so will always be available
        input_string = 'accnos={}'.format(self._tool_inputs['TSV_Accnos'][0])
        input_parameters = {'FASTA': 'fasta=',
                            'FASTQ': 'fastq=',
                            'TSV_Names': 'name=',
                            'TSV_Counts': 'count=',
                            'TSV_Groups': 'group=',
                            'TSV_AlignReport': 'alignreport=',
                            'TSV_List': 'list=',
                            'TSV_Taxonomy': 'taxonomy=',
                            'TSV_Qfile': 'qfile='}
        for key, input_files in self._tool_inputs.iteritems():
            # Only two keys are possible so the one that is not TSV_Accnos
            # will define the option flag that is needed
            if key != 'TSV_Accnos':
                input_string += ', ' + input_parameters[key]
                input_string += input_files[0].path
        input_string += ', outputdir={}'.format(self._folder)
        return input_string

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        # The extension of the output files is not always built in the same way. A dictionary
        # holds the suffixes that can be used when getting the base name
        output_extensions = {'FASTA': ['.', '.pick.fasta'],
                             'FASTQ': ['.', '.pick.fastq'],
                             'TSV_Names': ['.', '.pick.names'],
                             'TSV_Counts': ['.', '.pick.count_table'],
                             'TSV_Groups': ['.', '.pick.groups'],
                             'TSV_AlignReport': ['.align.report', '.pick.align.report'],
                             'TSV_List': ['.', '.pick.list'],
                             'TSV_Taxonomy': ['.', '.pick.taxonomy'],
                             'TSV_Qfile': ['.', '.pick.qual']}
        for key, input_files in self._tool_inputs.iteritems():
            # The TSV_Accnos file does not directly lead to output
            if key != 'TSV_Accnos':
                basename = super(MothurRemoveSeqs, self)._get_basename(key, output_extensions[key][0])
                self._tool_outputs[key] = [ToolIOFile(basename + output_extensions[key][1])]
