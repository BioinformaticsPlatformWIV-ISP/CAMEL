from pathlib import Path

SNAKEFILE_SCCMEC_TYPING = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_scc_mec_typing = Path('sccmec_typing')
OUTPUT_SCCMEC_TYPING_REPORT = _dir_scc_mec_typing / 'report' / 'html.io'
OUTPUT_SCCMEC_TYPING_REPORT_EMPTY = _dir_scc_mec_typing / 'report' / 'html-empty.io'
OUTPUT_SCCMEC_TYPING_SUMMARY = _dir_scc_mec_typing / 'summary_sccmec.tsv'
