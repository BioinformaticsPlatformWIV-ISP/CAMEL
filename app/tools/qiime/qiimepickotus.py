import os.path

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliofile import ToolIOFile
from app.tools.qiime.qiime import Qiime


class QiimePickOtus(Qiime):
    """
    The OTU picking step assigns similar sequences to operational taxonomic units, or OTUs, by clustering sequences
    based on a user-defined similarity threshold. Sequences which are similar at or above the threshold level are
    taken to represent the presence of a taxonomic unit (e.g., a genus, when the similarity threshold is set at 0.94)
    in the sequence collection.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(QiimePickOtus, self).__init__('qiime_pick_otus', '1.9.1', camel)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - FASTA_REF, BLAST_DB, SORTMERNA_DB are allowed as additional key
        - Only one file allowed per key
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Invalid input files (keys) given for pick_otus: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.iteritems():
            if key not in ['FASTA', 'FASTA_Ref', 'BLAST_DB', 'SORTMERNA_DB']:
                raise InvalidInputSpecificationError('Invalid input key given for pick_otus: {!r}'.format(self._tool_inputs))
            if len(self._tool_inputs[key]) != 1:
                raise InvalidInputSpecificationError('Invalid number (max = 1) of files in each key given for \
                                                     pick_otus: {!r}'.format(self._tool_inputs))

    def _set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        basename = super(QiimePickOtus, self)._get_basename()
        self._tool_outputs['TSV_Otu'] = [ToolIOFile(os.path.join(self._folder, basename + '_otus.txt'))]
        self._tool_outputs['LOG'] = [ToolIOFile(os.path.join(self._folder, basename + '_otus.log'))]

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        input_string = '-i {}'.format(self._tool_inputs['FASTA'][0])
        if os.path.isfile(os.path.join(self._folder, self._parameter_file)):
            input_string += ' -p {}'.format(os.path.join(self._folder, self._parameter_file))
        if 'FASTA_Ref' in self._tool_inputs:
            input_string += ' -r {}'.format(self._tool_inputs['FASTA_Ref'][0])
        if 'BLAST_DB' in self._tool_inputs:
            input_string += ' -b {}'.format(self._tool_inputs['BLAST_DB'][0])
        if 'SORTMERNA_DB' in self._tool_inputs:
            input_string += ' --sortmerna_db {}'.format(self._tool_inputs['SORTMERNA_DB'][0])
        input_string += ' -o {}'.format(self._folder)
        return input_string
