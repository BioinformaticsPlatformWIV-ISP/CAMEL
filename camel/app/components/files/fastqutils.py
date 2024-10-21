import gzip
import re
from itertools import zip_longest
from pathlib import Path
from typing import Set, Iterable, AbstractSet, BinaryIO, Dict, Union, Sequence

import screed
from Bio import SeqIO

from camel.app.command.command import Command
from camel.app.components.filesystemhelper import FileSystemHelper


class FastqUtils(object):

    """
    This class contains utility functions to work with FASTQ files.
    """

    @staticmethod
    def count_reads(infile: Path) -> int:
        """
        Count how many reads in a fastq file
        :param infile: file name of the fastq file to count
        :return: number of reads in fastq file
        """
        cat = 'zcat' if FileSystemHelper.is_gzipped(infile) else 'cat'
        cmd = f"{cat} {infile} | paste - - - - | wc -l"
        command = Command()
        command.command = cmd
        command.run(infile.resolve().parent)
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)
        return int(command.stdout.rstrip())

    @staticmethod
    def sort_fastq_by_identifier(infile: Path, outfile: Path, gzip_output: bool = False) -> None:
        """
        Function to sort the reads in a fastq file
        :param infile: file name of the file to sort
        :param outfile: file name of the sorted file
        :param gzip_output: gzip the output file
        :return: None
        """
        cat = 'zcat' if FileSystemHelper.is_gzipped(infile) else 'cat'
        gzip_ = '| gzip ' if gzip_output else ''
        cmd = f"{cat} {infile} | paste - - - - | sort -k1,1 -t \" \" | tr \"\t\" \"\n\" {gzip_}> {outfile}"
        command = Command()
        command.command = cmd
        command.run(infile.resolve().parent)
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)

    @staticmethod
    def convert_fastqs_to_interleaved_fastq(pe_1_file: Path, pe_2_file: Path, interleaved_file: Path,
                                            gzip_output: bool = False) -> None:
        """
        Convert two ordered fastq file into a interleaved fastq file, NO CHECKS are performed. If checks are needed,
        call 'create_paired_end' instead.
        :param pe_1_file: fastq file containing the first group of reads
        :param pe_2_file: fastq file containing the second group of reads
        :param interleaved_file: interleaved fastq file to be generated
        :param gzip_output: gzip the output file
        :return: None
        """
        # (z)cat is used as 'paste - - - - <(zcat filename)' does not work
        cat = 'zcat' if FileSystemHelper.is_gzipped(pe_1_file) else 'cat'
        gzip_ = '| gzip ' if gzip_output else ''
        cmd = f"paste <({cat} {pe_1_file} | paste - - - -) <({cat} {pe_2_file} | paste - - - -) | " \
              f"tr '\\t' '\\n' {gzip_}> {interleaved_file}"
        command = Command()
        command.command = cmd
        command.run(interleaved_file.resolve().parent)
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)

    @staticmethod
    def convert_interleaved_fastq_to_individual_fastqs(interleaved_file: Path, pe_1_file: Path, pe_2_file: Path,
                                                       gzip_output: bool = False) -> None:
        """
        Convert an interleaved fastq file into two fastq files each containing one group of reads. Input interleaved
        fastq file must be sorted. No CHECKS are performed.
        :param interleaved_file: interleaved fastq file containing both groups of reads
        :param pe_1_file: fastq file will hold the first group of reads
        :param pe_2_file: fastq file will hold the second group of reads
        :param gzip_output: gzip the output file
        :return: None
        """
        FastqUtils.split_interleaved_fastq(interleaved_file, pe_1_file, pe_2_file, gzip_output)

    @staticmethod
    def split_interleaved_fastq(interleaved_file: Path, pe_1_file: Path, pe_2_file: Path,
                                gzip_output: bool = False, pigz: bool = False) -> None:
        """
        Split an interleaved fastq file into two fastq files each containing one group of reads. Input interleaved fastq
        file must be sorted. No CHECKS are performed.
        :param interleaved_file: interleaved fastq file containing both groups of reads
        :param pe_1_file: fastq file will hold the first group of reads
        :param pe_2_file: fastq file will hold the second group of reads
        :param gzip_output: gzip the output file
        :param pigz: use pigz instead of gzip
        :return: None
        """
        cat = 'zcat' if FileSystemHelper.is_gzipped(interleaved_file) else 'cat'
        gzip_ = '| gzip ' if gzip_output and not pigz else '| pigz ' if gzip_output and pigz else ''
        cmd = f"{cat} {interleaved_file} | paste - - - - - - - - | tee >(cut -f 1-4 | tr '\t' '\n' {gzip_}> {pe_1_file}) |" \
              f" cut -f 5-8 | tr '\t' '\n' {gzip_}> {pe_2_file}"
        command = Command()
        command.command = cmd
        command.run(pe_1_file.resolve().parent)
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)

    @staticmethod
    def _get_read_name(read: screed.Record) -> str:
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
    def _write_read_to_file(read: screed.Record, outfile: BinaryIO):
        """
        Writes the given read to the given output file
        :param read: Read to write
        :param outfile: File to write to
        :return: None
        """
        outfile.write(f'@{read.name}\n{read.sequence}\n+\n{read.quality}\n'.encode(encoding='utf-8'))

    @staticmethod
    def _flush_se_reads(read_dict: Dict[str, screed.Record], outfile: BinaryIO):
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
    def _found_read(read: screed.Record, other_dict: Dict[str, screed.Record]) -> bool:
        """
        Checks whether the given read was already found in the dictionary of reads from the other file
        :param read: Read to check
        :param other_dict: Dictionary with the reads of the other file
        :return: True/False
        """
        return read is not None and FastqUtils._get_read_name(read) in other_dict

    @staticmethod
    def _process_reads(read: screed.Record, read_dict: Dict[str, screed.Record], other_dict: Dict[str, screed.Record],
                       pe_outf: BinaryIO, se_outf: BinaryIO) -> None:
        """
        Writes the paired reads to the paired output file and flushes the dictionary the read came from. As the files
        are sorted, the reads that came before the current read will be orphans and can be written to the SE file.
        :param read: Read to write
        :param read_dict: Dictionary with reads of the file the read came from
        :param other_dict: Dictionary with reads of the other file
        :param pe_outf: Output file handle for the paired reads
        :param se_outf: Output file handle for the orphans
        :return: None
        """
        FastqUtils._write_read_to_file(read_dict[FastqUtils._get_read_name(read)], pe_outf)
        read_dict.pop(FastqUtils._get_read_name(read))
        FastqUtils._write_read_to_file(other_dict[FastqUtils._get_read_name(read)], pe_outf)
        other_dict.pop(FastqUtils._get_read_name(read))
        FastqUtils._flush_se_reads(read_dict, se_outf)

    @staticmethod
    def create_paired_end(s1_file: Path, s2_file: Path, pe_out: Path, se_out: Path, gzip_output: bool = False) -> None:
        """
        Function to extract paired reads from two fastq files, i.e. fastq files with forward and reverse reads but
        of unequal length or non-paired reads. NOTE: The fastq files have to be sorted! Not sorting the files may result
        in reads being considered orphaned while they are not.
        :param s1_file: name of the first fastq file (forward reads - /1)
        :param s2_file: name of the second fastq file (reverse reads - /2)
        :param pe_out: file name for the paired end file (interleaved)
        :param se_out: file name for the orphan reads
        :param gzip_output: gzip the output file
        :return: None
        """
        read1_dict = {}
        read2_dict = {}
        outf_open_fn = gzip.open if gzip_output else open
        with outf_open_fn(pe_out, 'wb') as pe_outf, outf_open_fn(se_out, 'wb') as se_outf, \
                screed.open(s1_file) as screed_iter_1, screed.open(s2_file) as screed_iter_2:
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

    PATTERN_FQ_PE = r'(.+?)(_S\d+)?(_L\d{3})?[_.]R?1P?(_\d+)?.(fastq|fq)(.gz)?'
    PATTERN_FQ_SE = r'(.+?)(_S\d+)?(_L\d{3})?(_\d+)?.(fastq|fq)(.gz)?'

    @staticmethod
    def get_sample_name(fastq_path: Union[Path, str], pattern: str = PATTERN_FQ_PE) -> str:
        """
        Returns the sample name based on the given reads.
        :param fastq_path: FASTQ path
        :param pattern: Regex to determine the sample name
        :return: Sample name
        """
        basename = FileSystemHelper.make_valid(Path(fastq_path).name)
        m = re.match(pattern, basename, re.IGNORECASE)
        if m:
            return m.group(1)
        raise ValueError(f"Cannot determine sample name from: {basename}")

    @staticmethod
    def get_all_read_names(fastq_path: Path) -> Set[str]:
        """
        Retrieves all read names from the given fastq file
        :param fastq_path: Path to the fastq file
        :return: Set with read names
        """
        read_names = set()
        open_fn = gzip.open if FileSystemHelper.is_gzipped(fastq_path) else open
        with open_fn(fastq_path, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fastq'):
                read_names.add(record.id)
        return read_names

    @staticmethod
    def process_paired_end(pe_1_files: Sequence[Path], pe_2_files: Sequence[Path], se_files: Sequence[Path],
                           pe_out_1: Path, pe_out_2: Path, se_out: Path, gzip_output: bool = False) -> None:
        """
        Function to extract paired end reads from multiple (filtered) paired-end input files and output all
        other and orphaned reads to a single output file.
        :param pe_1_files: Iterable with forward read paired end files
        :param pe_2_files: Iterable with reverse read paired end files
        :param se_files: Iterable with single end files
        :param pe_out_1: File with all paired forward reads
        :param pe_out_2: File with all paired reverse reads
        :param se_out: File with all single end and orphaned reads
        :param gzip_output: gzip the output file
        :return: None
        """
        open_fn_in = gzip.open if FileSystemHelper.is_gzipped(pe_1_files[0]) else open
        open_fn_out = gzip.open if gzip_output else open
        fwd_reads = set()
        for file in pe_1_files:
            fwd_reads = fwd_reads | FastqUtils.get_all_read_names(file)
        rev_reads = set()
        for file in pe_2_files:
            rev_reads = rev_reads | FastqUtils.get_all_read_names(file)
        paired_reads = fwd_reads & rev_reads
        open_fn_out(se_out, 'wt').close()  # Initialize SE file so that append can be used without issues
        for infiles, outfile in [(pe_1_files, pe_out_1), (pe_2_files, pe_out_2)]:
            with open_fn_out(outfile, 'wt') as pe_outhandle, open_fn_out(se_out, 'at') as se_outhandle:
                for file in infiles:
                    for record in SeqIO.parse(open_fn_in(file, 'rt'), 'fastq'):
                        if record.id in paired_reads:
                            SeqIO.write(record, pe_outhandle, 'fastq')
                        else:
                            SeqIO.write(record, se_outhandle, 'fastq')
        with open_fn_out(se_out, 'at') as outhandle:
            for file in se_files:
                with open_fn_in(file, 'rt') as inhandle:
                    for line in inhandle:
                        outhandle.write(line)

    @staticmethod
    def process_paired_end_se(pe_1_files: Sequence[Path], pe_2_files: Sequence[Path], se_1_files: Sequence[Path], se_2_files: Sequence[Path],
                              pe_out_1: Path, pe_out_2: Path, se_out_1: Path, se_out_2: Path, gzip_output: bool = False) -> None:
        """
        Function to extract paired end reads from multiple (filtered) paired-end input files and output other and orphaned reads
        to single end output files.
        :param pe_1_files: Iterable with forward read paired end files
        :param pe_2_files: Iterable with reverse read paired end files
        :param se_1_files: Iterable with forward single end files
        :param se_2_files: Iterable with forward single end files
        :param pe_out_1: File with all paired forward reads
        :param pe_out_2: File with all paired reverse reads
        :param se_out_1: File with forward single end and orphaned reads
        :param se_out_2: File with reverse single end and orphaned reads
        :param gzip_output: gzip the output file
        :return: None
        """
        open_fn_in = gzip.open if FileSystemHelper.is_gzipped(pe_1_files[0]) else open
        open_fn_out = gzip.open if gzip_output else open
        fwd_reads = FastqUtils._get_read_names(pe_1_files)
        rev_reads = FastqUtils._get_read_names(pe_2_files)
        paired_reads = fwd_reads & rev_reads
        open_fn_out(se_out_1, 'w').close()  # Initialize SE file so that append can be used without issues
        open_fn_out(se_out_2, 'w').close()
        for infiles, outfile_pe, outfile_se in [(pe_1_files, pe_out_1, se_out_1), (pe_2_files, pe_out_2, se_out_2)]:
            with open_fn_out(outfile_pe, 'wt') as pe_outhandle, open_fn_out(outfile_se, 'at') as se_outhandle:
                for file in infiles:
                    for record in SeqIO.parse(open_fn_in(file, 'rt'), 'fastq'):
                        if record.id in paired_reads:
                            SeqIO.write(record, pe_outhandle, 'fastq')
                        else:
                            SeqIO.write(record, se_outhandle, 'fastq')
        for infile, outfile in [(se_1_files, se_out_1), (se_2_files, se_out_2)]:
            with open_fn_out(outfile, 'at') as outhandle:
                for file in infile:
                    with open_fn_in(file, 'rt') as inhandle:
                        for line in inhandle:
                            outhandle.write(line)

    @staticmethod
    def _get_read_names(files: Iterable[Path]) -> AbstractSet[str]:
        """
        Returns a set with all the read names in the given files
        :param files: Files for which the read names need to be extracted
        :return: Set of read names
        """
        read_names = set()
        for file in files:
            read_names = read_names | FastqUtils.get_all_read_names(file)
        return read_names

    @staticmethod
    def count_bases(input_file: Path) -> int:
        """
        Calculates the number of bases in the given input files
        :param input_file: File path
        :return: Number of bases
        """
        cat = 'zcat' if FileSystemHelper.is_gzipped(input_file) else 'cat'
        cmd = f"{cat} {input_file} | paste - - - - | cut -f 2 | tr -d '\n' | wc -c"
        command = Command()
        command.command = cmd
        command.run(Path.cwd())
        return int(command.stdout)
