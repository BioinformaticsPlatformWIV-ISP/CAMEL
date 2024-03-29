from pathlib import Path

SNAKEFILE_PREPROCESS = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_pre_process = Path('preprocess')

OUTPUT_PRE_PROCESS_REPORT = _dir_pre_process / 'report' / 'html.io'
OUTPUT_PRE_PROCESS_SUMMARY = _dir_pre_process / 'report' / 'summary_preprocess.tsv'
OUTPUT_PRE_PROCESS_INFORMS = _dir_pre_process / 'report' / 'informs.io'
OUTPUT_PRE_PROCESS_FASTQ = _dir_pre_process / 'downsample' / 'fq_dict.io'

# AmpliGone
OUTPUT_PRE_PROCESS_AMPLIGONE_REPORT = _dir_pre_process / 'ampligone' / 'html.io'
OUTPUT_PRE_PROCESS_AMPLIGONE_REPORT_EMPTY = _dir_pre_process / 'ampligone' / 'html-empty.io'

# Amplicon clipping
OUTPUT_PRE_PROCESS_CLIPPING_REPORT = _dir_pre_process / 'ampliconclip' / 'html.io'
OUTPUT_PRE_PROCESS_CLIPPING_REPORT_EMPTY = _dir_pre_process / 'ampliconclip' / 'html-empty.io'
