import os
import re

from app.error.invalidparametererror import InvalidParameterError
from app.io.tooliofile import ToolIOFile
from app.tools.bowtie2.bowtie2 import Bowtie2


class Bowtie2Map(Bowtie2):
    """
    Reads mapping using Bowtie2. It does **not** support using both PE and SE reads. Does **not** support interleaved
    fastq format due to lack of use.
    """
    OUTPUT_NAME = 'bowtie2_readmap.sam'
    SAMPLE_NAME = 'sampleA'

    def __init__(self, camel):
        """
        Initialize Bowtie2
        :param camel: Camel instance
        :return: None
        """
        super(Bowtie2Map, self).__init__('bowtie2 map', '2.3.0', camel)

        self._time_inform_mapping = {
            'Time loading reference': 'time_loading_reference',
            'Time loading forward index': 'time_loading_forward_index',
            'Time loading mirror index': 'time_loading_reverse_index',
            'Time searching': 'time_searching',
            'Multiseed full-index search': 'time_multiseed_full-index_search',
            'Overall time': 'time_total'
        }
        self._align_inform_mapping = {
            'were unpaired; ': 'stats_singe_reads_in',
            'were paired; ': 'stats_paired_reads_in',
            'aligned concordantly 0 times': 'stats_pair_0_concord_map',
            'aligned concordantly exactly 1 time': 'stats_pair_1_concord_map',
            'aligned concordantly >1 times': 'stats_pair_n_concord_maps',
            'aligned discordantly 1 time': 'stats_pair_disconcord_maps',
            'aligned 0 times': 'stats_single_0_map',
            'aligned exactly 1 time': 'stats_single_1_map',
            'aligned >1 times': 'stats_single_n_maps',
            'reads; of these:': 'stats_reads_count',
            'overall alignment rate': 'stats_map_rate'
        }
        self._mod = None
        self._fastq_inputs_str = ''
        self._refgenome_str = ''
        self._output_str = ''
        self._readgroup_str = ''

    def _execute_tool(self):
        """
        Function to run Bowtie2 to map reads
        :return: None
        """
        self.__set_input()
        self.__set_output()
        self.__build_command()
        self._execute_command()
        self.__set_inform()

    def __set_input(self):
        """
        Set extra input
        :return: None
        """
        if 'SAMPLE_NAME' in self._tool_inputs:
            sample_name = self._tool_inputs['SAMPLE_NAME'][0].value
        else:
            sample_name = Bowtie2Map.SAMPLE_NAME
        self._readgroup_str += " --rg-id {!r}".format(sample_name)

    def _check_input(self):
        """
        Check input for Bowtie2 mapping
        :return: None
        """
        # Note that Bowtie2 can map both PE and SE reads together
        if 'FASTQ_PE' in self._tool_inputs:
            self._mod = 'PE'
            if len(self._tool_inputs['FASTQ_PE']) != 2:
                raise ValueError("Paired end input requires exactly 2 files.")
            self._fastq_inputs_str += ' -1 {} -2 {}'.format(
                self._tool_inputs['FASTQ_PE'][0].path, self._tool_inputs['FASTQ_PE'][1].path
            )
        if 'FASTQ_SE' in self._tool_inputs:
            if not self._mod:
                self._mod = 'SE'
            self._fastq_inputs_str += ' -U {}'.format(",".join(f.path for f in self._tool_inputs['FASTQ_SE']))
        if not self._fastq_inputs_str:
            raise ValueError("No FASTQ_PE of FASTQ_SE input found")

        if 'INDEX_GENOME_PREFIX' not in self._tool_inputs:
            raise ValueError('No genome index input (INDEX_GENOME_PREFIX) found.')
        self._refgenome_str = "-x {}".format(self._tool_inputs['INDEX_GENOME_PREFIX'][0].value)

        super(Bowtie2Map, self)._check_input()

    def __set_output(self):
        """
        Set output for Bowtie2 read mapping
        :return None
        """
        sam_filename = os.path.join(self._folder, Bowtie2Map.OUTPUT_NAME)
        self._tool_outputs['SAM'] = [ToolIOFile(sam_filename)]
        self._output_str = "-S {}".format(sam_filename)

    @staticmethod
    def __check_mode_exclusiveness(options):
        """
        Alignment mode exclusiveness check
        :param options: names of commandline options set
        """
        if ('end-to-end' in options) and ('local' in options):
            raise InvalidParameterError(
                "Bowtie2 reads mapping modes 'end-to-end' and 'local' are exclusive to each other, cannot be specifed at the same time!"
            )

    @staticmethod
    def __check_presets_exclusiveness(align_mode, presets, options):
        """
        Check whether more than one preset is specified for a given alignment mode
        :param align_mode: current alignment mode
        :param presets: presets of current alignment mode
        :param options: names of commandline options set
        :return: None
        """
        p_count = 0
        for p in presets:
            if p in options:
                p_count += 1
        if p_count > 1:
            raise InvalidParameterError(
                "Cannot set more than one preset for Bowtie2 mode {!r}, parameters set: {}".format(
                    align_mode, ",".join(options))
            )

    @staticmethod
    def __check_wrong_preset(align_mode, options, wrong_presets):
        """
        Check whether wrong preset is specified for a given alignment mode
        :param align_mode: current alignment mode
        :param options: names of commandline options set
        :param wrong_presets: presets should not be used with current alignment mode
        :return: None
        """
        for p in wrong_presets:
            if p in options:
                raise InvalidParameterError(
                    "Bowtie2 incompatible preset: reads mapping mode {}, preset {!r}.".format(align_mode, p))

    @staticmethod
    def __check_mode_preset_conflicts(options):
        """
        Alignment mode and preset conflicts check
        :param options: names of set options
        """
        local_preset = [
            'very-fast-local',
            'fast-local',
            'sensitive-local',
            'very-sensitive-local'
        ]
        end_to_end_preset = [
            'very-fast',
            'fast',
            'sensitive',
            'very-sensitive'
        ]
        if 'end-to-end' in options:
            Bowtie2Map.__check_wrong_preset('end-to-end', options, local_preset)
            Bowtie2Map.__check_presets_exclusiveness('end-to-end', end_to_end_preset, options)
        elif 'local' in options:
            Bowtie2Map.__check_wrong_preset('local', options, end_to_end_preset)
            Bowtie2Map.__check_presets_exclusiveness('local', local_preset, options)

    def _check_parameters(self):
        """
        Check the exclusiveness of different mods of Bowtie2 and the exclusiveness of presets for each mod and between mods
        """
        super(Bowtie2Map, self)._check_parameters()

        options = set(self._parameters.keys())
        Bowtie2Map.__check_mode_exclusiveness(options)
        Bowtie2Map.__check_mode_preset_conflicts(options)

    def __build_command(self):
        """
        Build command to run Bowtie2
        :return: None
        """
        self._command.command = '{} {} {} {} {} {}'.format(
            self._tool_command,
            " ".join(self._build_options()),
            self._readgroup_str,
            self._refgenome_str,
            self._fastq_inputs_str,
            self._output_str
        )

    def __set_time_info(self, line):
        """
        Set running time related information
        :param line: the content of current line (from self._command.stderr)
        :return: boolean, True if a time information is found
        """
        # time format: starts with number and composed of number and ":"
        res = re.search(': (\d[\d|:]+)', line)
        if res:
            time = res.group().replace(": ", "")
            for pattern, key in self._time_inform_mapping.items():
                if line.find(pattern) >= 0:
                    self.informs[key] = time
                    return True

        return False

    def __set_mapping_inform(self, line):
        """
        Set mapping result statistics information
        :param line: the content of current line (from self._command.stderr)
        :return: None
        """
        # search for all numbers
        res = re.findall("([\d|.]+)", line)
        if res:
            if line.find("pairs aligned 0 times concordantly or discordantly") > 0:
                # NOTE: this one also matches "aligned 0 times", but this is more specific, hence
                #       should be checked first
                self.informs['stats_pair_map_single'] = "{}".format(res[0])
                return

            for pattern, key in self._align_inform_mapping.items():
                if line.find(pattern) > 0:
                    if len(res) > 1:
                        # two numbers found: count, percentage
                        self.informs[key] = "{} ({}%)".format(res[0], res[1])
                    else:
                        # only one number found, count/percentage
                        self.informs[key] = "{}".format(res[0])
                    return

    def __set_inform(self):
        """
        Analyse the result of Bowtie2 reads mapping, and extra result statistics into tool inform
        :return: None
        """
        self.informs['tool_name'] = 'Bowtie2'
        self.informs['mod'] = self._mod

        # parse output to extract information
        for l in self.stderr.splitlines():
            time_inform_set = self.__set_time_info(l)
            if not time_inform_set:
                self.__set_mapping_inform(l)
