from pathlib import Path

SNAKEFILE_GENOMETYPING = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_genometyping = Path('genometyping')
OUTPUT_GENOMETYPING_REPORT = _dir_genometyping / 'report' / 'html.io'
OUTPUT_GENOMETYPING_SUMMARY = _dir_genometyping / 'summary_out.tsv'
OUTPUT_GENOMETYPING_FASTA_REF = _dir_genometyping / 'blastn_processing' / 'fasta.io'
OUTPUT_GENOMETYPING_INDEX_GENOME_PREFIX = _dir_genometyping / 'blastn_processing' / 'index_genome_prefix.io'
OUTPUT_SEQTK_SUBSAMPLE_FASTQ = _dir_genometyping / 'seqtk_subsample' / 'fastq.io'
OUTPUT_SEQTK_SUBSAMPLE_INFORMS = _dir_genometyping / 'seqtk_subsample' / 'informs.io'
OUTPUT_SEQTK_CONVERT_FASTA = _dir_genometyping / 'seqtk_convert' / 'fasta.io'
OUTPUT_SEQTK_CONVERT_INFORMS = _dir_genometyping / 'seqtk_convert' / 'informs.io'
OUTPUT_BLASTN_ASN = _dir_genometyping / 'blastn' / 'asn.io'
OUTPUT_BLASTN_TSV = _dir_genometyping / 'blastn' / 'tsv.io'
OUTPUT_BLASTN_INFORMS = _dir_genometyping / 'blastn' / 'informs.io'
OUTPUT_BLASTN_PROCESSING_INFORMS = _dir_genometyping / 'blastn_processing' / 'informs.io'
