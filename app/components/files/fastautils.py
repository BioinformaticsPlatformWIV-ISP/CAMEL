import os
from Bio import SeqIO
from app.command.command import Command


class FastaUtils(object):

    """
    Helper to perform FASTA file related functions
    """

    @staticmethod
    def read_as_index_dict(fasta):
        """
        Read in fasta file as an index dictionary of SeqRecord keyed by sequence id. Handle big files, as the complete
        sequence is retrieved when access a SeqRecord.
        :param fasta: fasta file to read
        :return: a sequence record dictionary
        """
        return SeqIO.index(fasta, "fasta")

    @staticmethod
    def read_as_dict(fasta):
        """
        Read in fasta file as a dictionary keyed by sequence id. More efficient for small files.
        :param fasta: fasta file to read
        :return: a sequence record dictionary
        """
        return SeqIO.to_dict(SeqIO.parse(open(fasta, 'rU'), "fasta"))

    @staticmethod
    def write(sequences, fasta):
        """
        Write a list of sequence records into fasta file
        :param sequences: list of sequences as SeqRecord object
        :param fasta: fasta file to write into
        :return: None
        """
        with open(fasta, 'w') as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")

    @staticmethod
    def count_reads(infile):
        """
        Count how many reads in a fastq file
        :param infile: file name of the fasta file to count
        :return: number of reads in fasta file
        """
        cmd = "cat {!r} | paste - - | wc -l".format(infile)
        command = Command()
        command.command = cmd
        command.run_command(os.path.dirname(os.path.abspath(infile)))
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)
        return int(command.stdout.rstrip())
