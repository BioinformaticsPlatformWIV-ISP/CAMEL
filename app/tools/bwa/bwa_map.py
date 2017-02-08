import os

from app.io.tooliofile import ToolIOFile
from app.tools.bwa.bwa import BWA


class BWAMap(BWA):
    """
    Reads mapping using 'bwa mem' from BWA

    UPDATE:
    - since v0.7.15, bwa mem support both PE and SE reads mapping through interleaved fastq file with revamped '-p'
      option (smart pairing) (see changelog)
    """
    OUTPUT_NAME = 'bwa_readmap.sam'
    DEFAULT_SAMPLE_NAME = 'sampleA'

    def __init__(self, camel):
        """
        Initialize BWAMap
        :param camel: Camel instance
        :return: None
        """
        super(BWAMap, self).__init__('bwa_mem', '0.7.15', camel)
        self._fastq_inputs_str = None
        self._readgroup_str = ''

    def _execute_tool(self):
        """
        Function to run BWA mem to map reads
        :return: None
        """
        self.__set_input()
        self.__set_output()
        self.__build_command()
        self._execute_command()

    def __set_input(self):
        """
        Set input
        :return: None
        """
        if 'FASTQ_PE' in self._tool_inputs:
            self._fastq_inputs_str = ' '.join(f.path for f in self._tool_inputs['FASTQ_PE'])
        elif 'FASTQ_SE' in self._tool_inputs:
            self._fastq_inputs_str = self._tool_inputs['FASTQ_SE'][0].path
        elif 'FASTQ_INT' in self._tool_inputs:
            self._fastq_inputs_str = "-p {}".format(self._tool_inputs['FASTQ_INT'][0].path)

        if 'SAMPLE_NAME' in self._tool_inputs:
            sample_name = self._tool_inputs['SAMPLE_NAME'][0].value
        else:
            sample_name = BWAMap.DEFAULT_SAMPLE_NAME
        # Read Group format: '@RG\tID:foo\tSM:bar'
        self._readgroup_str += "@RG\tID:{0}\tSM:{0}".format(sample_name)

    def _check_input(self):
        """
        Check input for BWA mem.
        :return: None
        """
        super(BWAMap, self)._check_input()

        if 'FASTQ_PE' in self._tool_inputs:
            if len(self._tool_inputs['FASTQ_PE']) != 2:
                raise ValueError("Paired end fastq inputs require exactly 2 files.")
        elif 'FASTQ_SE' in self._tool_inputs:
            if len(self._tool_inputs['FASTQ_SE']) > 1:
                raise ValueError("Single end fastq input requires exactly 1 file.")
        elif 'FASTQ_INT' in self._tool_inputs:
            if len(self._tool_inputs['FASTQ_INT']) > 1:
                raise ValueError("Interleaved fastq input requires exactly 1 file.")
        else:
            raise ValueError("No supported FASTQ_PE, FASTQ_SE, or FASTQ_INT input found")

        if 'INDEX_GENOME_PREFIX' not in self._tool_inputs:
            raise ValueError('No genome index input (INDEX_GENOME_PREFIX) found.')

    def __set_output(self):
        """
        Set proper outputs for BAW mem
        :return: None
        """
        self._tool_outputs['SAM'] = [ToolIOFile(os.path.join(self._folder, BWAMap.OUTPUT_NAME))]

    def __build_command(self):
        """
        Build command to run BWA mem
        :return: None
        """
        self._command.command = '{} {} -R {!r} {} {} > {}'.format(
            self._tool_command,
            ' '.join(self._build_options()),
            self._readgroup_str,
            self._tool_inputs['INDEX_GENOME_PREFIX'][0].value,
            self._fastq_inputs_str,
            self._tool_outputs['SAM'][0].path
        )
