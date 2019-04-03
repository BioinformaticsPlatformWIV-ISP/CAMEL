import os.path

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurPcrSeqs(Mothur):
    """
    The pcr.seqs will trim inputted sequences based on a variety of user-defined options.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('mothur_pcr_seqs', '1.39.1', camel)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA is required
        - Only TSV_Oligos is allowed as additional key
        - Only one FASTA file is allowed
        :return: None
        """
        super(MothurPcrSeqs, self)._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('No input file given for Mothur pcr.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidInputSpecificationError('Invalid number (max = 1) of files given for Mothur \
                                                 pcr.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) > 2:
            raise InvalidInputSpecificationError('Too many input keys given for Mothur pcr.seqs: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'TSV_Oligos']:
                raise InvalidInputSpecificationError('Invalid input key given for Mothur pcr.seqs: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = ['fasta={}'.format(self._tool_inputs['FASTA'][0])]
        if 'TSV_Oligos' in self._tool_inputs:
            items.append('oligos={}'.format(self._tool_inputs['TSV_Oligos'][0]))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = super(MothurPcrSeqs, self)._get_basename()
        extension = self._tool_inputs['FASTA'][0].file_extension
        self._tool_outputs['FASTA'] = [ToolIOFile('{}.pcr{}'.format(basename, extension))]
        if os.path.isfile('{}.bad.accnos'.format(basename)):
            self._tool_outputs['TEXT'] = [ToolIOFile('{}.bad.accnos'.format(basename))]

    def _check_command_output(self):
        """
        Analyzes output to discover if the run was successful. If an error was present in stdout, a RuntimeError is
        raised and stdout is displayed
        :return: None
        """
        for line in self.stdout.splitlines():
            if line.startswith('[ERROR]: name mismatch in pcr.seqs'):
                # Hopefully temporary fix for bug in Mothur that gives these error messages
                pass
            elif line.startswith('[ERROR]') or line.startswith('Unable to open'):
                raise RuntimeError(self.stdout + '\n' + '!!! Mothur failed to run !!! See above for more information.')
