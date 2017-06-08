import os
import pysam


class SAMBAMutils(object):

    """
    Helper to perform SAM BAM file related functions
    """

    @staticmethod
    def get_record_count(infile):
        """
        Count the total number of records in a SAM/BAM file
        :param infile: input SAM/BAM file
        :return: the number of records in SAM/BAM file
        """
        _, file_ext = os.path.splitext(infile)
        if file_ext.lower() == '.sam':
            map_file = pysam.AlignmentFile(infile, "r")
        elif file_ext.lower() == '.bam':
            map_file = pysam.AlignmentFile(infile, "rb")

        return map_file.count(until_eof=True)

    @staticmethod
    def is_empty(infile):
        """
        Check whether a SAM/BAM file contains only header
        :param infile: input SAM/BAM file
        :return: True if infile contains no read records
        """
        return SAMBAMutils.get_record_count(infile) == 0
