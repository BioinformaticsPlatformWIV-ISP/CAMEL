from pathlib import Path

SNAKEFILE_QC = f'{Path(__file__).parent / Path(__file__).stem}.smk'
SNAKEFILE_QC_summary = f'{Path(__file__).parent / Path(__file__).stem}_summary.smk'
