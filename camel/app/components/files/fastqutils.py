from itertools import zip_longest

import os
import re
import screed


from camel.app.command.command import Command


class FastqUtils(object):

    """
    Helper to perform FASTQ file related functions
    """
    # Note that func convert_interleaved_fastq_to_individual_fastqs has been renamed as split_interleaved_fastq

    def __init__(self):
        pass

    @staticmethod
    def count_reads(infile):
        """
        Count how many reads in a fastq file
        :param infile: file name of the fastq file to count
        :return: number of reads in fastq file
        """
        cmd = "cat {!r} | paste - - - - | wc -l".format(infile)
        command = Command()
        command.command = cmd
        command.run_command(os.path.dirname(os.path.abspath(infile)))
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)
        return int(command.stdout.rstrip())

    @staticmethod
    def sort_fastq_by_identifier(infile, outfile):
        """
        Function to sort the reads in a fastq file
        :param infile: file name of the file to sort
        :param outfile: file name of the sorted file
        :return: None
        """
        cmd = "cat {!r} | paste - - - - | sort -k1,1 -t \" \" | tr \"\t\" \"\n\" > {!r}".format(infile, outfile)
        command = Command()
        command.command = cmd
        command.run_command(os.path.dirname(os.path.abspath(outfile)))
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)

    @staticmethod
    def convert_fastqs_to_interleaved_fastq(pe_1_file, pe_2_file, interleaved_file):
        """
        Convert two ordered fastq file into a interleaved fastq file, NO CHECKS are performed. If checks are needed, call
        'create_paired_end' instead.

        :param pe_1_file: fastq file containing the first group of reads
        :param pe_2_file: fastq file containing the second group of reads
        :param interleaved_file: interleaved fastq file to be generated
        :return: None
        """
        cmd = "paste <(paste - - - - < {!r}) <(paste - - - - < {!r}) | tr '\t' '\n' \ > {!r}".format(
            pe_1_file, pe_2_file, interleaved_file)
        command = Command()
        command.command = cmd
        command.run_command(os.path.dirname(os.path.abspath(interleaved_file)))
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)

    @staticmethod
    def convert_interleaved_fastq_to_individual_fastqs(interleaved_file, pe_1_file, pe_2_file):
        """
        Convert a interleaved fastq file into two fastq files each containing one group of reads. Input interleaved fastq
        file must be sorted. No CHECKS are performed.
        :param interleaved_file: interleaved fastq file containing both groups of reads
        :param pe_1_file: fastq file will hold the first group of reads
        :param pe_2_file: fastq file will hold the second group of reads
        :return: None
        """
        FastqUtils.split_interleaved_fastq(interleaved_file, pe_1_file, pe_2_file)

    @staticmethod
    def split_interleaved_fastq(interleaved_file, pe_1_file, pe_2_file):
        """
        Split a interleaved fastq file into two fastq files each containing one group of reads. Input interleaved fastq
        file must be sorted. No CHECKS are performed.
        :param interleaved_file: interleaved fastq file containing both groups of reads
        :param pe_1_file: fastq file will hold the first group of reads
        :param pe_2_file: fastq file will hold the second group of reads
        :return: None
        """
        cmd = "paste - - - - - - - - < {!r} | tee >(cut -f 1-4 | tr '\t' '\n' > {!r}) | cut -f 5-8 | tr '\t' '\n' > {!r}".format(
            interleaved_file, pe_1_file, pe_2_file)
        command = Command()
        command.command = cmd
        command.run_command(os.path.dirname(os.path.abspath(pe_1_file)))
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)

    @staticmethod
    def _get_read_name(read):
        """
        Splits the read identification at the first whitespace character. Handles the cases: 'name/1' - '@name 1:rst' -
        '@name seq/1'.
        :param read: The read to process
        :return: Read name
        """
        parts = read.name.split(None, 1)
        lhs, rhs = [parts[0], parts[1] if len(parts) > 1 else '']
        if lhs.endswith('/1') or lhs.endswith('/2'):
            return lhs.split('/', 1)[0]
        elif rhs.startswith('1:') or rhs.startswith('2:'):
            return lhs
        elif rhs.endswith('/1') or rhs.endswith('/2'):
            return rhs.split('/', 1)[0]
        else:
            raise Exception("Read name is in an illegal format: {}".format(read.name))

    @staticmethod
    def _write_read_to_file(read_out, outfile):
        """
        Writes the given read to the given output file
        :param read_out: Read to write
        :param outfile: File to write to
        :return: None
        """
        outstr = '@{name}\n{sequence}\n+\n{quality}\n'.format(
            name=read_out.name,
            sequence=read_out.sequence,
            quality=read_out.quality)
        outfile.write(outstr)

    @staticmethod
    def _flush_se_reads(read_dict, outfile):
        """
        Remove all the reads from a dictionary after writing them to a file. As the files are sorted these will be the
        orphaned reads.
        :param read_dict: Dictionary containing reads
        :param outfile: File to write to
        :return:
        """
        for key in list(read_dict.keys()):
            FastqUtils._write_read_to_file(read_dict[key], outfile)
            read_dict.pop(key)

    @staticmethod
    def _found_read(read, other_dict):
        """
        Checks whether the given read was already found in the dictionary of reads from the other file
        :param read: Read to check
        :param other_dict: Dictionary with the reads of the other file
        :return: True/False
        """
        if read is not None and FastqUtils._get_read_name(read) in other_dict:
            return True
        return False

    @staticmethod
    def _process_reads(read, read_dict, other_dict, pe_outf, se_outf):
        """
        Writes the paired reads to the paired output file and flushes the dictionary the read came from. As the files
        are sorted, the reads that came before the current read will be orphans and can be written to the SE file.
        :param read: Read to write
        :param read_dict: Dictionary with reads of the file the read came from
        :param other_dict: Dictionary with reads of the other file
        :param pe_outf: Output file for the paired reads
        :param se_outf: Output file for the orphans
        :return: None
        """
        FastqUtils._write_read_to_file(read_dict[FastqUtils._get_read_name(read)], pe_outf)
        read_dict.pop(FastqUtils._get_read_name(read))
        FastqUtils._write_read_to_file(other_dict[FastqUtils._get_read_name(read)], pe_outf)
        other_dict.pop(FastqUtils._get_read_name(read))
        FastqUtils._flush_se_reads(read_dict, se_outf)

    @staticmethod
    def create_paired_end(s1_file, s2_file, pe_out, se_out):
        """
        Function to extract paired reads from two fastq files, i.e. fastq files with forward and reverse reads but
        of unequal length or non-paired reads. NOTE: The fastq files have to be sorted! Not sorting the files may result
        in reads being considered orphaned while they are not.
        :param s1_file: name of the first fastq file (forward reads - /1)
        :param s2_file: name of the second fastq file (reverse reads - /2)
        :param pe_out: file name for the paired end file (interleaved)
        :param se_out: file name for the orphan reads
        """
        read1_dict = {}
        read2_dict = {}
        with open(pe_out, 'wb') as pe_outf, open(se_out, 'wb') as se_outf, screed.open(s1_file) as screed_iter_1, screed.open(s2_file) as screed_iter_2:
            for read1, read2 in zip_longest(screed_iter_1, screed_iter_2):
                # When the end of one file is reached before the other one, the line from that file will be None
                if read1 is not None:
                    read1_dict[FastqUtils._get_read_name(read1)] = read1
                if read2 is not None:
                    read2_dict[FastqUtils._get_read_name(read2)] = read2

                if FastqUtils._found_read(read1, read2_dict):
                    FastqUtils._process_reads(read1, read1_dict, read2_dict, pe_outf, se_outf)

                if FastqUtils._found_read(read2, read1_dict):
                    FastqUtils._process_reads(read2, read2_dict, read1_dict, pe_outf, se_outf)

            FastqUtils._flush_se_reads(read1_dict, se_outf)
            FastqUtils._flush_se_reads(read2_dict, se_outf)

    @staticmethod
    def get_sample_name(fastq_filename):
        """
        Returns the sample name based on the given reads. It tries to match the following formats(in this order):
        - Sample - Name_S\d + _L\d + _R[12]_\d + .fastq(e.g.: S15BD00757_S20_L001_R2_001.fastq)
        - Sample - Name_1.fastq, Sample - Name_2.fastq(e.g.: reads_1.fastq)
        :
            param fastq_filename:
                FASTQ filename
        :
            return:
                Sample name
        """
        m = re.match('(.*)_S\d+_L\d+_R[12]_\d+.[fastq]+$', fastq_filename)
        if m:
            return m.group(1)
        m = re.match('(.*)_[12].[fastq]+$', fastq_filename)
        if m:
            return m.group(1)
        raise ValueError("Cannot determine sample name from: {}".format(fastq_filename))
