from pathlib import Path

from camel.resources.snakefile import assembly_flye

SNAKEFILE_QC = f'{Path(__file__).parent / Path(__file__).stem}.smk'

consensus_by_tool = {
    'Flye':  assembly_flye.OUTPUT_ASSEMBLY_FASTA,
    'Medaka': Path('polish', 'long_reads', 'flye', 'sequence', 'fasta.io'),
    'Polypolish': Path('polish', 'short_reads', 'flye', 'polypolish', 'fasta.io'),
    'POLCA': Path('polish', 'short_reads', 'flye', 'polca', 'fasta.io'),
    'Unicycler': Path('unicycler', 'fasta.io'),
    'Medaka-Unicycler': Path('polish', 'long_reads', 'unicycler', 'sequence', 'fasta.io'),
    'Polypolish-Unicycler': Path('polish', 'short_reads', 'unicycler', 'polypolish', 'fasta.io'),
    'POLCA-Unicycler': Path('polish', 'short_reads', 'unicycler', 'polca', 'fasta.io'),
}

ALE_KEYS = ['depth', 'kmer', 'insert', 'place']
