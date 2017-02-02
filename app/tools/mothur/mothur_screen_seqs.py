from app.tools.mothur.mothur import Mothur
from app.io.tooliofile import ToolIOFile


class MothurScreenSeqs(Mothur):
    """
    The screen.seqs command enables you to keep sequences that fulfill certain user defined criteria. Furthermore, it
    enables you to cull those sequences not meeting the criteria from a names, group, contigsreport, alignreport and
    summary file.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(MothurScreenSeqs, self).__init__('mothur_screen_seqs', '1.39.0', camel)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - Only 'FASTA', 'TSV_Summary', 'TSV_Groups', 'TSV_AlignReport', 'TSV_ContigReport',
          'TSV_Names', 'TSV_Counts', 'TSV_Qfile', and 'TSV_Taxonomy' are allowed
        - Only one input file per key is allowed
        :return: None
        """
        super(MothurScreenSeqs, self)._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise RuntimeError('No input file given for Mothur screen.seqs: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.iteritems():
            if key not in ['FASTA', 'TSV_Summary', 'TSV_Groups', 'TSV_AlignReport', 'TSV_ContigReport',
                           'TSV_Names', 'TSV_Counts', 'TSV_Qfile', 'TSV_Taxonomy']:
                raise RuntimeError('Invalid input key given for Mothur screen.seqs: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise RuntimeError('Invalid number (max = 1) of files given for Mothur \
                                   screen.seqs: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and input/output directories
        Example: ffastq=FileR1.fastq, rfastq=FileR2.fastq, inputdir=/test/data/input/,
        outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        # Start with empty string as the keys may not be ordered in self._tool_inputs
        input_string = ''
        # As there can be only one file per key, this first file of the list is added
        for key, input_files in self._tool_inputs.iteritems():
            if key == 'FASTA':
                input_string += 'fasta={}, '.format(input_files[0])
            elif key == 'TSV_Groups':
                input_string += 'group={}, '.format(input_files[0])
            elif key == 'TSV_Summary':
                input_string += 'summary={}, '.format(input_files[0])
            elif key == 'TSV_Names':
                input_string += 'name={}, '.format(input_files[0])
            elif key == 'TSV_AlignReport':
                input_string += 'alignreport={}, '.format(input_files[0])
            elif key == 'TSV_ContigsReport':
                input_string += 'contigreport={}, '.format(input_files[0])
            elif key == 'TSV_Taxonomy':
                input_string += 'taxonomy={}, '.format(input_files[0])
            elif key == 'TSV_Counts':
                input_string += 'count={}, '.format(input_files[0])
        # As the input string will now end with ', ' this will have to be removed
        input_string = input_string[:-2]
        input_string += ', outputdir={}'.format(self._folder)
        return input_string

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        # Screen.seqs re-appends the original extension (e.g. .fasta or .align) after first adding .good
        path = self._tool_inputs['FASTA'][0].path
        fasta_extension = path[self._tool_inputs['FASTA'][0].path.rfind('.'):]
        self._tool_outputs['FASTA'] = [ToolIOFile(super(MothurScreenSeqs, self)._get_basename('FASTA') + '.good' + fasta_extension)]
        self._tool_outputs['TSV_Bad'] = [ToolIOFile(super(MothurScreenSeqs, self)._get_basename('FASTA') + '.bad.accnos')]
        # Depending on the input keys more files will be created
        if 'TSV_Groups' in self._tool_inputs:
            self._tool_outputs['TSV_Groups'] = [ToolIOFile(super(MothurScreenSeqs, self)._get_basename('TSV_Groups', '.groups') +
                                                '.good.groups')]
        if 'TSV_Summary' in self._tool_inputs:
            self._tool_outputs['TSV_Summary'] = [ToolIOFile(super(MothurScreenSeqs, self)._get_basename('TSV_Summary', '.summary') +
                                                 '.good.summary')]
        if 'TSV_Names' in self._tool_inputs:
            self._tool_outputs['TSV_Names'] = [ToolIOFile(super(MothurScreenSeqs, self)._get_basename('TSV_Names', '.names') +
                                               '.good.names')]
        if 'TSV_AlignReport' in self._tool_inputs:
            self._tool_outputs['TSV_AlignReport'] = [ToolIOFile(super(MothurScreenSeqs, self)._get_basename('TSV_AlignReport', '.align.report') +
                                                     '.good.align.report')]
        if 'TSV_ContigsReport' in self._tool_inputs:
            self._tool_outputs['TSV_ContigsReport'] = [ToolIOFile(super(MothurScreenSeqs, self)._get_basename('TSV_ContigReport', '.contig.report') +
                                                       '.good.contigs.report')]
        if 'TSV_Counts' in self._tool_inputs:
            self._tool_outputs['TSV_Counts'] = [ToolIOFile(super(MothurScreenSeqs, self)._get_basename('TSV_Counts', '.count_table') +
                                                '.good.count_table')]
