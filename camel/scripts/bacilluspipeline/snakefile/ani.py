from pathlib import Path

SNAKEFILE_ANI = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_ani = Path('ani')
INPUT_FASTA_ANI = _dir_ani / 'fasta.io'
OUTPUT_VAL_ANI = _dir_ani / 'val-ani.io'
OUTPUT_INFORMS_ANI = _dir_ani / 'informs.io'
OUTPUT_ANI_REPORT = _dir_ani / 'html.io'
OUTPUT_ANI_REPORT_EMPTY = _dir_ani / 'html-empty.io'
OUTPUT_ANI_SUMMARY = _dir_ani / 'summary_out.tsv'
