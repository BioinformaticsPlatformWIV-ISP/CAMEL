from pathlib import Path

SNAKEFILE_SUBTYPING = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_subtyping = Path('subtyping')
OUTPUT_SUBTYPING_REPORT = _dir_subtyping / 'report' / 'html.io'
OUTPUT_SUBTYPING_SUMMARY = _dir_subtyping / 'summary' / 'summary_out.tsv'
OUTPUT_SEQTK_SUBSAMPLE_FASTQ = _dir_subtyping / 'seqtk_subsample' / 'fastq.io'
OUTPUT_SEQTK_SUBSAMPLE_INFORMS = _dir_subtyping / 'seqtk_subsample' / 'informs.io'
OUTPUT_SEQTK_CONVERT_FASTA = _dir_subtyping / 'seqtk_convert' / 'fasta.io'
OUTPUT_SEQTK_CONVERT_INFORMS = _dir_subtyping / 'seqtk_convert' / 'informs.io'
OUTPUT_BLASTN_ASN = _dir_subtyping / 'blastn' / 'asn.io'
OUTPUT_BLASTN_INFORMS = _dir_subtyping / 'blastn' / 'informs.io'
OUTPUT_BLASTN_PROCESSING_INFORMS = _dir_subtyping / 'blastn_processing' / 'informs.io'
