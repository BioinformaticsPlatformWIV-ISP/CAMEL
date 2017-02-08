import re


class FastqUtils(object):
    """
    Class with utility functions for FASTQ files.
    """

    @staticmethod
    def get_sample_name(fastq_filename):
        """
        Returns the sample name based on the given reads. It tries to match the following formats (in this order):
        - Sample-Name_S\d+_L\d+_R[12]_\d+.fastq (e.g.: S15BD00757_S20_L001_R2_001.fastq)
        - Sample-Name_1.fastq, Sample-Name_2.fastq (e.g.: reads_1.fastq)
        :param fastq_filename: FASTQ filename
        :return: Sample name
        """
        m = re.match('(.*)_S\d+_L\d+_R[12]_\d+.[fastq]+$', fastq_filename)
        if m:
            return m.group(1)
        m = re.match('(.*)_[12].[fastq]+$', fastq_filename)
        if m:
            return m.group(1)
        raise ValueError("Cannot determine sample name from: {}".format(fastq_filename))
