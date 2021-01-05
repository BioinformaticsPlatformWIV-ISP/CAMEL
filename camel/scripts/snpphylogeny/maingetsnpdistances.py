import pandas as pd
import argparse
from typing import Optional, Sequence

from Bio import SeqIO


def parse_args(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """
    Parses the command line arguments.
    :param args: Arguments to parse
    :return: Parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--output', required=True, help='Output file with SNP distances')
    return parser.parse_args(args)


def calc_snp_distance(seq_a: str, seq_b: str) -> int:
    """
    Calculates the SNP distances between two sequences.
    N's are ignored in the calculation.
    :param seq_a: Sequence A
    :param seq_b: Sequence B
    :return: None
    """
    distance = 0
    for nuc_a, nuc_b in zip(seq_a, seq_b):
        if nuc_a == 'N' or nuc_b == 'N':
            continue
        distance += 1 if nuc_a != nuc_b else 0
    return distance


if __name__ == '__main__':
    options = parse_args()

    # Parse input
    with open(options.fasta) as handle:
        seqs = SeqIO.parse(handle, 'fasta')
        seq_by_id = {s.id: str(s.seq) for s in seqs}

    # Calculate distances
    data_distance = pd.DataFrame([{'sample': id_} for id_ in sorted(seq_by_id.keys())])
    for sample in sorted(seq_by_id.keys()):
        data_distance[sample] = data_distance.apply(
            lambda x: calc_snp_distance(seq_by_id[x['sample']], seq_by_id[sample]), axis=1)

    # Save output
    data_distance.to_csv(options.output, sep='\t', index=False)
