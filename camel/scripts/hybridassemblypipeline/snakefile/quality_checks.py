from pathlib import Path

SNAKEFILE_QC = f'{Path(__file__).parent / Path(__file__).stem}.smk'
consensus_by_tool = {
    'Flye': Path('assembly_flye', 'filtering', 'assembly_filtered.fasta'),
    'Medaka': Path('medaka', 'consensus.fasta'),
    'Polypolish': Path('polishing', 'polypolish', 'polished.fasta'),
    'POLCA': Path('polishing', 'polca', 'polished.fasta'),
    'Unicycler': Path('unicycler', 'assembly.fasta')
}

ALE_KEYS = ['depth', 'kmer', 'insert', 'place']
