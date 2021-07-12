from pathlib import Path

SNAKEFILE_ALIGNMENT = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_alignment = Path('alignment')
OUTPUT_ALIGNMENT_BAM_UNMAPPED = _dir_alignment/ "fastq_to_ubam" / "{input_basename}.unmapped.bam.io"
OUTPUT_ALIGNMENT_BAM = _dir_alignment / "gather_bqsr_sorted_bam" / "bqsr_gathered_sorted.bam.io"
OUTPUT_MARK_DUPLICATES_METRICS = _dir_alignment / 'mark_duplicates' / "duplicate_metrics.txt.io"