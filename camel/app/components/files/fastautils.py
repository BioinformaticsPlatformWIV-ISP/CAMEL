import math
from pathlib import Path
from typing import List, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from camel.app.command.command import Command


class FastaUtils(object):

    """
    Helper to perform FASTA file related functions
    """

    @staticmethod
    def read_as_index_dict(fasta: Path):
        """
        Read in fasta file as an index dictionary of SeqRecord keyed by sequence id. Handle big files, as the complete
        sequence is retrieved when access a SeqRecord.
        :param fasta: fasta file to read
        :return: a sequence record dictionary
        """
        return SeqIO.index(fasta, "fasta")

    @staticmethod
    def read_as_dict(fasta: Path) -> dict:
        """
        Read in fasta file as a dictionary keyed by sequence id. More efficient for small files.
        :param fasta: fasta file to read
        :return: a sequence record dictionary
        """
        with fasta.open() as handle:
            return SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    @staticmethod
    def write(sequences: List[Union[SeqIO.SeqRecord, str]], fasta: Path) -> None:
        """
        Write a list of sequence records into fasta file
        :param sequences: list of sequences as SeqRecord object
        :param fasta: fasta file to write into
        :return: None
        """
        with open(fasta, 'w') as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")

    @staticmethod
    def count_reads(infile: Path) -> int:
        """
        Count how many reads in a fasta file
        :param infile: file name of the fasta file to count
        :return: number of reads in fasta file
        """
        cmd = f'grep -c "^>" {infile}'
        command = Command()
        command.command = cmd
        command.run(infile.resolve().parent)
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)
        return int(command.stdout.rstrip())

    @staticmethod
    def batch_iterator(iterator, batch_size: int) -> List[any]:
        """
        Returns lists of length batch_size. This is a generator function, and it
        returns lists of the entries from the supplied iterator.  Each list will
        have batch_size entries, although the final list may be shorter.
        :param iterator: Iterator to create batches from
        :param batch_size: Size of each batch
        :return: List from iterator with batch_size entries
        """
        entry = True  # Make sure we loop once
        while entry:
            batch = []
            while len(batch) < batch_size:
                try:
                    entry = iterator.next()
                except StopIteration:
                    entry = None
                if entry is None:
                    # End of file
                    break
                batch.append(entry)
            if batch:
                yield batch

    @staticmethod
    def split_fasta(fasta: Path, outdir: Path, n_parts: int = None, parts_size: int = None) -> List[str]:
        """
        Splits a fasta file in the given number of parts or in parts of the given size.
        :param fasta: Fasta file to split
        :param outdir: Output directory where the split files can be placed
        :param n_parts: Number of parts to split into
        :param parts_size: Size of the parts that need to be created
        :return: List with the filenames of the parts that are created
        """
        if n_parts and parts_size:
            raise ValueError('Number of parts and parts size are mutually exclusive!')
        if n_parts:
            n_reads = FastaUtils.count_reads(fasta)
            parts_size = math.ceil(n_reads/n_parts)
        basename = fasta.stem
        groups = []
        record_iter = SeqIO.parse(open(fasta), 'fasta')
        for i, batch in enumerate(FastaUtils.batch_iterator(record_iter, parts_size)):
            filename = f'{outdir / basename}.group{i+1}.fasta'
            with open(filename, 'w') as outhandle:
                SeqIO.write(batch, outhandle, 'fasta')
            groups.append(filename)
        return groups

    @staticmethod
    def is_indexed(fasta: Path) -> bool:
        """
        Checks if the input fasta file is indexed with the FAI index.
        :param fasta: Fasta file to check for associated index
        :return: True or False depending on the presence of the index file
        """
        try:
            next(fasta.parent.glob(f'{fasta.name}.fai'))
        except StopIteration:
            return False
        return True

    @staticmethod
    def convert_fasta_to_fastq(fasta_input_file: Path, fastq_output_file: Path) -> None:
        """
        Converts a FASTA file to a FASTQ file with perfect Phred quality scores.
        :return: None
        """
        with fasta_input_file.open('r') as fasta_file, fastq_output_file.open('w') as fastq_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                # Create a fake quality score string of the same length as the sequence
                fake_quality = [40] * len(record.seq)

                # Create a SeqRecord with the same sequence and a fake quality string
                fake_record = SeqRecord(record.seq, id=record.id, description=record.description, letter_annotations={
                    "phred_quality": fake_quality})
                # Write the SeqRecord in FASTQ format
                SeqIO.write(fake_record, fastq_file, 'fastq')
