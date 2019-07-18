from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REPEATS_FASTA = '/data/pipelines/saureus/sparepeats.fasta'
PROFILES_TSV = '/data/pipelines/saureus/spatypes.csv'

if __name__ == '__main__':
    with open(REPEATS_FASTA) as handle:
        seqs = list(SeqIO.parse(handle, 'fasta'))
        seq_by_number = {s.id: s for s in seqs}

    profiles = {}
    with open(PROFILES_TSV) as handle:
        for line in handle.readlines():
            if line.startswith('NT'):
                continue
            profile, repeats = line.strip().split(',')
            profiles[profile] = [int(x) for x in repeats.split('-')]

    profiles_to_check = ['t002', 't011', 't018', 't128', 't437', 't529']
    # profiles_to_check = ['t529']
    for profile in profiles_to_check:
        repeat_nbs = profiles[profile]
        sequence = ''.join(str(seq_by_number[f'r{n:02d}'].seq) for n in repeat_nbs)

        print(f'>{profile}')
        print(sequence)

        print(f'>{profile}_rc')
        print(str(Seq(sequence).reverse_complement()))

        print(f'>{profile}_rc-step')
        print(''.join(str(seq_by_number[f'r{n:02d}'].seq.reverse_complement()) for n in reversed(repeat_nbs)))

    # Create FASTA file with all sequences
    all_seqs = []
    for profile, repeats in sorted(profiles.items()):
        s = Seq(''.join(str(seq_by_number[f'r{n:02d}'].seq) for n in repeats))
        all_seqs.append(SeqRecord(s, profile, description=''))
    with open('/scratch/bebog/working/profiles.fasta', 'w') as handle:
        SeqIO.write(all_seqs, handle, 'fasta')

