from pathlib import Path

SNAKEFILE_MAPPING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_align = Path('alignment')
OUTPUT_ALIGNMENT_SAM = _dir_align / 'sam.io'
OUTPUT_ALIGNMENT_INFORMS = _dir_align / 'informs.io'
OUTPUT_ALIGNMENT_SUMMARY = _dir_align / 'summary.tsv'
OUTPUT_ALIGNMENT_BAM = _dir_align / 'bam.io'
OUTPUT_ALIGNMENT_ALIGNMENTSUMMARY = _dir_align / 'alignmentsummary.io'
OUTPUT_ALIGNMENT_INSERTSIZE = _dir_align / 'txt_insertsize.io'
OUTPUT_ALIGNMENT_MAPQUALITYDISTRIBUTION = _dir_align / 'txt_mapqualitydistribution.io'
OUTPUT_ALIGNMENT_MAPQUALITYDISTRIBUTION_PDF = _dir_align / 'pdf_mapqualitydistribution.io'
OUTPUT_ALIGNMENT_GCBIAS = _dir_align / 'gcbias.io'
OUTPUT_ALIGNMENT_GCBIAS_SUMMARY = _dir_align / 'gcbias_summary.io'
OUTPUT_ALIGNMENT_GCBIAS_FIGURE = _dir_align / 'gcbias_figure.io'
OUTPUT_ALIGNMENT_SAMTOOLS_DEPTH = _dir_align / 'samtools_depth_tsv.io'
OUTPUT_ALIGNMENT_SAMTOOLS_DEPTH_INFORMS = _dir_align / 'samtools_depth_informs.io'
OUTPUT_ALIGNMENT_SAMTOOLS_DEPTH_ANALYZER_INFORMS = _dir_align / 'samtools_depth_analyzer_informs.io'
OUTPUT_ALIGNMENT_PICARD_METRICS = _dir_align / 'informs_picard_metrics.io'
OUTPUT_ALIGNMENT_REPORT = _dir_align / 'report' / 'html.io'

