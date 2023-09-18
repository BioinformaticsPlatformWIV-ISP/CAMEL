from pathlib import Path

SNAKEFILE_QC = f'{Path(__file__).parent / Path(__file__).stem}.smk'

consensus_by_tool = {
    'Flye': Path('assembly_flye', 'filtering', 'fasta.io'),
    'Medaka': Path('medaka', 'flye', 'stitch', 'fasta.io'),
    'Polypolish': Path('polishing', 'flye', 'polypolish', 'polished.fasta'),
    'POLCA': Path('polishing', 'flye', 'polca', 'fasta.io'),
    'Unicycler': Path('unicycler', 'assembly.fasta'),
    'Medaka-Unicycler': Path('medaka', 'unicycler', 'stitch', 'fasta.io'),
    'Polypolish-Unicycler': Path('polishing', 'unicycler', 'polypolish', 'polished.fasta'),
    'POLCA-Unicycler': Path('polishing', 'unicycler', 'polca', 'fasta.io'),
}

ALE_KEYS = ['depth', 'kmer', 'insert', 'place']
