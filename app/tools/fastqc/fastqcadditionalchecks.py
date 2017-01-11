import logging

from app.tools.tool import Tool


class FastQCAdditionalChecks(Tool):
    """
    Class that performs various quality checks based on FastQC output.
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: Camel instance
        """
        super(FastQCAdditionalChecks, self).__init__('FastQC additional checks', '0.1', camel)

    def _check_input(self):
        """
        Checks whether all the required input files are in the inputs.
        :return: None
        """
        if 'TXT' not in self._tool_inputs:
            raise ValueError("No TXT input found.")
        super(FastQCAdditionalChecks, self)._check_input()

    def _execute_tool(self):
        """
        Performs the quality checks.
        :return: None
        """
        for input_file in self._tool_inputs['TXT']:
            modules = self.__get_modules(input_file)
            current_informs = {
                'Mean Q-score drop':
                    self.__test_mean_qscore_drop(modules['Per base sequence quality']),

                'Average quality score':
                    self.__test_average_read_quality(modules['Per sequence quality scores']),

                'Per base sequence content':
                    self.__test_per_base_sequence_content(modules['Per base sequence content']),

                'GC content':
                    self.__test_gc_content(modules['Basic Statistics']),

                'Maximal N-fraction':
                    self.__test_max_n_count(modules['Per base N content']),

                'Sequence length distribution':
                    self.__test_sequence_length_distribution(modules['Sequence Length Distribution'])
            }
            self.informs[self.__get_sample_name(modules['Basic Statistics'])] = current_informs

    @staticmethod
    def __get_modules(input_file):
        """
        Get the content of the different modules of the fastqc data.
        :param input_file: Name of the fastqc data file
        :return: Modules with the content
        """
        modules = {}
        with open(input_file.path) as file_:
            key = None
            content = []
            for line in file_.readlines():
                if line.strip() == '>>END_MODULE':
                    modules[key] = content
                    content = []
                elif line.startswith('>>'):
                    key = line[2:].strip().split('\t')[0]
                else:
                    content.append(line.strip())
        return modules

    @staticmethod
    def __get_sample_name(data):
        """
        Returns the sample name.
        :param data: Basic statistics data
        :return: GC content
        """
        for line in data:
            if line.startswith('Filename'):
                return line.split('\t')[1]

    @staticmethod
    def __get_mean_qscore_drop(data, threshold):
        """
        Returns the base where the mean qscore drops below the threshold.
        :param data: Per base sequence quality data
        :param threshold: Threshold
        :return: base
        """
        for row in data[1:]:
            base = row.split('\t')[0]
            mean_q_score = row.split('\t')[1]
            if float(mean_q_score) < threshold:
                return float(base.split('-')[0])
        return float('inf')

    def __test_mean_qscore_drop(self, data):
        """
        Test to see at which base the mean qscore drops below a threshold.
        :param data: Mean qscore data
        :return: 'Pass', 'Warning' or 'Fail'.
        """
        fail_length = float(self._parameters['qscore_drop_fail_length'].value)
        warn_length = float(self._parameters['qscore_drop_warn_length'].value)
        threshold = float(self._parameters['qscore_drop_threshold'].value)

        qscore_drop = self.__get_mean_qscore_drop(data, threshold)
        logging.debug("Mean qscore drop at {}".format(qscore_drop))

        if qscore_drop < fail_length:
            return 'Fail'
        elif qscore_drop < warn_length:
            return 'Warn'
        else:
            return 'Pass'

    @staticmethod
    def __get_average_read_quality(data):
        """
        Returns the average read quality.
        :param data: Per sequence quality scores data
        :return: Average read quality
        """
        total_count = 0.0
        total_quality = 0.0
        for row in data[1:]:
            quality, count = row.split('\t')
            total_count += float(count)
            total_quality += float(quality) * float(count)
        return total_quality / total_count

    def __test_average_read_quality(self, data):
        """
        Tests the average read quality.
        :param data: Per sequence quality scores data
        :return: 'Pass', 'Warning' or 'Fail'
        """
        read_quality_fail = float(self._parameters['average_read_quality_fail'].value)
        read_quality_warn = float(self._parameters['average_read_quality_warn'].value)
        average_read_quality = self.__get_average_read_quality(data)
        logging.debug("Average read quality: {:.2f}".format(average_read_quality))

        if average_read_quality < read_quality_fail:
            return 'Fail'
        elif average_read_quality < read_quality_warn:
            return 'Warn'
        else:
            return 'Pass'

    @staticmethod
    def __get_max_sequence_content_difference(data, nb_of_skipped_bases, nb_of_skipped_bases_end):
        """
        Returns the maximal difference between A-T and C-G.
        :param data: Per base sequence content data
        :param nb_of_skipped_bases: Nb of bases to skip in the start of the reads
        :param nb_of_skipped_bases: Nb of bases to skip at the end of the reads
        :return: Difference
        """
        max_difference = 0.0
        last_base = int(data[-1].split('\t')[0].split('-')[1])
        logging.debug("Checking A-T, G-C difference between base {} and {}".format(
            nb_of_skipped_bases, last_base - nb_of_skipped_bases_end))
        for row in data[1:]:
            base, freq_g, freq_a, freq_t, freq_c = row.split('\t')
            interval_upper = int(base.split('-')[-1])
            if nb_of_skipped_bases < interval_upper <= (last_base - nb_of_skipped_bases_end):
                at_difference = abs(float(freq_a) - float(freq_t))
                if at_difference > max_difference:
                    max_difference = at_difference
                gc_difference = abs(float(freq_c) - float(freq_g))
                if gc_difference > max_difference:
                    max_difference = gc_difference
        return max_difference

    def __test_per_base_sequence_content(self, data):
        """
        Test whether difference between A-T & C-G is below a threshold at every position.
        :param data: Per base sequence content data
        :return: 'Pass', 'Warn' or 'Fail'
        """
        per_base_sequence_content_fail = float(self._parameters['per_base_sequence_content_fail'].value)
        per_base_sequence_content_warn = float(self._parameters['per_base_sequence_content_warn'].value)
        skipped_bases = int(self._parameters['per_base_sequence_content_skipped'].value)
        skipped_bases_end = int(self._parameters['per_base_sequence_content_skipped_end'].value)

        max_difference = self.__get_max_sequence_content_difference(data, skipped_bases, skipped_bases_end)
        logging.debug("Maximal difference between A/T or C/G: {:.2f}".format(max_difference))
        if max_difference > per_base_sequence_content_fail:
            return 'Fail'
        elif max_difference > per_base_sequence_content_warn:
            return 'Warn'
        else:
            return 'Pass'

    @staticmethod
    def __get_gc_content_median(data):
        """
        Returns the median of the GC content.
        :param data: Per sequence GC content data
        :return: Median GC content
        """
        total_count = 0.0
        for row in data[1:]:
            _, count = row.split('\t')
            total_count += float(count)
        middle = total_count / 2

        sum_count = 0.0
        prev_gc = 0
        for row in data[1:]:
            gc, count = row.split('\t')
            sum_count += float(count)
            if sum_count > float(middle):
                return int(prev_gc)
            prev_gc = gc

    @staticmethod
    def __get_gc_content_modus(data):
        """
        Returns the modus of the GC content.
        :param data: Per sequence GC content data
        :return: Modus GC content
        """
        max_count = 0.0
        modus_gc = 0
        for row in data[1:]:
            gc, count = row.split('\t')
            if float(count) > max_count:
                max_count = float(count)
                modus_gc = int(gc)
        return modus_gc

    @staticmethod
    def __get_gc_content(data):
        """
        Returns the GC content.
        :param data: Basic statistics data
        :return: GC content
        """
        for line in data:
            if line.startswith('%GC'):
                return float(line.split('\t')[1])

    def __test_gc_content(self, data):
        """
        Checks the GC content modus.
        :param data: Per sequence GC content data
        :return: 'Pass', 'Warn' or 'Fail'
        """
        reference_gc_content = float(self._parameters['gc_content_reference'].value)
        gc_content_difference_warn = float(self._parameters['gc_content_difference_warn'].value)
        gc_content_difference_fail = float(self._parameters['gc_content_difference_fail'].value)

        gc_content = self.__get_gc_content(data)
        logging.debug("Detected GC content: {:.2f}".format(gc_content))
        difference = abs(gc_content - reference_gc_content)

        if difference > gc_content_difference_fail:
            return 'Fail'
        elif difference > gc_content_difference_warn:
            return 'Warn'
        else:
            return 'Pass'

    @staticmethod
    def __get_max_n_count(data):
        """
        Returns the maximum N count.
        :param data: Per base N content data
        :return: Count
        """
        max_count = 0.0
        for row in data[1:]:
            base, n_count = row.split('\t')
            if float(n_count) > max_count:
                max_count = float(n_count)
        return max_count

    def __test_max_n_count(self, data):
        """
        Tests whether the N count is below a threshold for every base.
        :param data: Per base N content data
        :return: 'Pass', 'Warn' or 'Fail'
        """
        threshold_fail = float(self._parameters['n_count_threshold_fail'].value)
        threshold_warn = float(self._parameters['n_count_threshold_warn'].value)
        max_n_count = self.__get_max_n_count(data)
        logging.debug("Maximal N count: {:.4f}".format(max_n_count))

        if max_n_count > threshold_fail:
            return 'Fail'
        elif max_n_count > threshold_warn:
            return 'Warn'
        else:
            return 'Pass'

    @staticmethod
    def __get_fraction_of_sequences_below_threshold(data, threshold):
        """
        Checks whether the given percentage of read lengths if above the given threshold.
        :param data: Sequence length distribution data
        :param threshold: Threshold
        :return: Percentage of reads below threshold
        """
        sequences_below_threshold = 0.0
        total_sequences = 0.0
        for row in data[1:]:
            length_interval, count = row.split('\t')
            interval_lower = float(length_interval.split('-')[1])
            if interval_lower < threshold:
                sequences_below_threshold += float(count)
            total_sequences += float(count)
        return sequences_below_threshold / total_sequences

    def __test_sequence_length_distribution(self, data):
        """
        Checks if there are not too much short sequences.
        :param data: Sequence length distribution data
        :return: 'Pass', 'Warn' or 'Fail'
        """
        threshold_fail = float(self._parameters['sequence_length_threshold_fail'].value)
        threshold_warn = float(self._parameters['sequence_length_threshold_warn'].value)
        max_fraction = float(self._parameters['sequence_length_fraction'].value)

        fraction_below_fail = self.__get_fraction_of_sequences_below_threshold(data, threshold_fail)
        logging.debug("Fraction of sequences with length below fail threshold: {:.2f}".format(fraction_below_fail))

        fraction_below_warn = self.__get_fraction_of_sequences_below_threshold(data, threshold_warn)
        logging.debug("Fraction of sequences with length below warn threshold: {:.2f}".format(fraction_below_warn))

        if self.__get_fraction_of_sequences_below_threshold(data, threshold_fail) > max_fraction:
            return 'Fail'
        elif self.__get_fraction_of_sequences_below_threshold(data, threshold_warn) > max_fraction:
            return 'Warn'
        else:
            return 'Pass'
