from pathlib import Path

SNAKEFILE_QC = f'{Path(__file__).parent / Path(__file__).stem}.smk'
consensus_by_tool = {
    'Flye': Path('assembly_flye', 'filtering', 'assembly_filtered.fasta'),
    'Medaka': Path('medaka', 'consensus.fasta'),
    'POLCA': Path('polishing', 'polca', 'polished.fasta'),
    'Polypolish': Path('polishing', 'polypolish', 'polished.fasta'),
    'Unicycler': Path('unicycler', 'assembly.fasta')
}