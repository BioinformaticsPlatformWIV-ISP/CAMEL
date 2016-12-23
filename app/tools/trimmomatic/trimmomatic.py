import os
import re

from app.tools.tool import Tool


class Trimmomatic(Tool):
    """
    A flexible read trimming tool for Illumina NGS Data.
    """

    def __init__(self, camel):
        """
        Initializes Trimmomatic.
        :return: None
        """
        super(Trimmomatic, self).__init__('Trimmomatic', '0.32', camel)
        self._mod = ''
        self._input_string = None
        self._output_string = None
        self.baseout = None

    def _check_output(self, folder):
        """
        Check the resulting output files and fills the common stream
        object with them - empty output file(s) are removed from the
        reported outputs, hence expect empty list
        :return: None
        """
        output_files = None
        if self._mod == 'PE':
            basename = os.path.splitext(self.baseout)[0]
            output_files = {
                'FASTQ_PE': [os.path.join(folder, basename + '_1P.fastq'),
                             os.path.join(folder, basename + '_2P.fastq')],
                'FASTQ_SE_FORWARD': [os.path.join(folder, basename + '_1U.fastq')],
                'FASTQ_SE_REVERSE': [os.path.join(folder, basename + '_2U.fastq')]
            }
            # remove empty output files
            if os.stat(output_files['FASTQ_PE'][0]).st_size == 0:
                # for paired reads output, if one file is empty, both are empty!
                output_files['FASTQ_PE'] = []
            if os.stat(output_files['FASTQ_SE_FORWARD'][0]).st_size == 0:
                output_files['FASTQ_SE_FORWARD'] = []
            if os.stat(output_files['FASTQ_SE_REVERSE'][0]).st_size == 0:
                output_files['FASTQ_SE_REVERSE'] = []
        elif self._mod == 'SE':
            if os.stat(self.baseout).st_size == 0:
                raise StandardError(
                    "No reads pass Trimmomatic QC check! Input: {}".format(
                        self._tool_inputs)
                )
            output_files = {'FASTQ': self.baseout}

        self._tool_outputs = output_files

    def _check_input(self):
        """
        Check inputs and set the mod of trimmomatic, only accept fastq files, two mod: 'pe' for paired-end reads, 'se' 
        for single-end reads
        :return: None
        """
        super(Trimmomatic, self)._check_input()
        self.baseout = self._parameters['baseout'].value
        if re.search('\.gz$', self.baseout):
            raise RuntimeError(
                'System require trimmomatic outputs plain unziped fastq file, a ziped name specified: {}'.format(
                    self.baseout)
            )
        if 'FASTQ_PE' in self._tool_inputs:
            if len(self._tool_inputs['FASTQ_PE']) != 2:
                raise RuntimeError('FASTQ_PE specified as inputs for trimmomatic, not exactly 2 files found.'
                                   ' tool_inputs[FASTQ_PE]: {!r}'.format(self._tool_inputs['FASTQ_PE']))
            self._mod = 'PE'
            self._input_string = " ".join([in_file.path for in_file in self._tool_inputs['FASTQ_PE']])
        elif 'FASTQ_SE' in self._tool_inputs:
            if len(self._tool_inputs['FASTQ_SE']) != 1:
                raise RuntimeError('FASTQ_SE specified as inputs for trimmomatic, not exactly 1 file found.'
                                   ' tool_inputs[FASTQ_SE]: {!r}'.format(self._tool_inputs['FASTQ_SE']))
            self._mod = 'SE'
            self._input_string = self._tool_inputs['FASTQ_SE'][0].path
        else:
            raise RuntimeError('Required input files does not found for trimmomatic, '
                               'require FASTQ_PE or FASTQ_SE found none.')

    def _build_pe_command(self, trim_opts):
        """
        Set PE mod command
        :param trim_opts: basic trimming options
        :return: None
        """
        self._output_string = '-baseout {}'.format(self.baseout)
        #
        # thread options
        thread_opt = str(self._parameters['threads'])

        # trimming option: illumina pe library prep adapter
        trim_opts.insert(0, self._parameters['illuminaclip-PE'].option + self._parameters['illuminaclip-PE'].value)
        trim_opt_string = " ".join(trim_opts)
        #
        # NOTE: when use -baseout option, output speicification
        #       must come BEFORE the input specified as individual
        #       files
        self._command.command = " ".join([
            self._tool_command, self._mod, thread_opt, self._output_string, self._input_string, trim_opt_string
        ])

    def _build_se_command(self, trim_opts):
        """
        Set SE mod command
        :param trim_opts: basic trimming options
        :return: None
        """
        self._output_string = self.baseout
        #
        # thread options
        thread_opt = " ".join(self._parameters['threads'])
        #
        # trimming option: illumina se library prep adapter
        trim_opts.insert(0, "".join(self._parameters['illuminaclip-SE']))
        trim_opt_string = " ".join(trim_opts)
        #
        # compose command
        self._command.command = " ".join([
            self._tool_command, self._mod, thread_opt, self._input_string, self._output_string, trim_opt_string
        ])

    def build_command(self):
        """
        Concatenates required parameters and options to build the command to run trimmomatic
        :return: None
        """
        # mod independent options
        #
        # trimming options
        excluded_parameters = ['baseout', 'threads',
                               'illuminaclip-PE', 'illuminaclip-SE']
        trim_opts = self._build_options(excluded_parameters, '')

        # mod dependent options
        if self._mod == 'PE':
            self._build_pe_command(trim_opts)

        elif self._mod == 'SE':
            self._build_se_command(trim_opts)

    def analyze_result(self):
        """
        Analyzes output to discover if run was successful and fill self.inform with result statistics
        :return: status of trimmomatic run (True: succeed/False: failed)
        """
        succeed = False
        try:
            informs = self.informs
            informs['mod'] = self._mod

            for line in self._command.stdout.splitlines():
                qc_encoding = re.search(
                    'Quality encoding detected as (?P<encode>\w+)', line)
                if qc_encoding:
                    # store input file encoding information in cfgs (which
                    # to be stored in the process ini file)
                    informs['encoding'] = qc_encoding.group('encode')
                elif re.match('Input Read Pairs', line):
                    # FASTQ_PE as input
                    # 'Input Read Pairs: 42000 Both Surviving: 41269 (98.26%) Forward Only Surviving: 729 (1.74%)
                    # Reverse Only Surviving: 0 (0.00%) Dropped: 2 (0.00%)'
                    res = re.match(r'Input Read Pairs: (\d+) Both Surviving: (\d+ \([\d\.%]+\)) '
                                   r'Forward Only Surviving: (\d+ \([\d\.%]+\)) '
                                   r'Reverse Only Surviving: (\d+ \([\d\.%]+\)) '
                                   r'Dropped: (\d+ \([\d\.%]+\))', line)
                    informs['paired_reads_in'] = res.groups()[0]
                    informs['paired_reads_out'] = res.groups()[1]
                    informs['forward_only_reads'] = res.groups()[2]
                    informs['reverse_only_reads'] = res.groups()[3]
                    informs['reads_drop'] = res.groups()[4]

                elif re.match('Input Read', line):
                    # FASTQ_SE as input
                    # 'Input Reads: 42000 Surviving: 41608 (99.07%) Dropped: 392 (0.93%)'
                    res = re.match(r'Input Reads: (\d+) '
                                   r'Surviving: (\d+ \([\d\.%]+\)) '
                                   r'Dropped: (\d+ \([\d\.%]+\))', line)
                    informs['reads_in'] = res.groups()[0]
                    informs['reads_out'] = res.groups()[1]
                    informs['reads_drop'] = res.groups()[2]

                elif line.strip() == 'Exit status: 0':
                    succeed = True

            if not succeed:
                raise RuntimeError(
                    'Error running trimmomatic. Please check the standard output and the stderr for more information.')

        except AttributeError:
            # self._logger.error(
            #     "%s. No splitlines attribute in output", self.name)
            print('No splitlines output')

    def run(self, folder='.'):
        """
        Function to run Trimmomatic tool
        :param folder: Folder where the command is executed.
        :return: None
        """
        self._check_input()
        self.build_command()
        super(Trimmomatic, self).run(folder)
        self.analyze_result()
        self._check_output(folder)
