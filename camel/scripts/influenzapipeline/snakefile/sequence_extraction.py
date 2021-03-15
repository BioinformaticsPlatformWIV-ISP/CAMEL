from pathlib import Path

SNAKEFILE_SEQ_EXTRACTION = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_seq_extraction = Path('seq_extraction')
OUTPUT_SEQ_EXTRACTION_ADD_OR_REPLACE_BAM = _dir_seq_extraction / 'add_or_replace' / 'bam.io'
OUTPUT_SEQ_EXTRACTION_FASTQ_TO_SAM = _dir_seq_extraction / 'fastq_to_sam' / 'bam.io'
OUTPUT_SEQ_EXTRACTION_MERGE_BAM_ALIGNMENT = _dir_seq_extraction / 'merge_bam_alignment' / 'bam.io'
OUTPUT_SEQ_EXTRACTION_MARK_DUPLICATES = _dir_seq_extraction / 'mark_duplicates' / 'bam.io'
OUTPUT_SEQ_EXTRACTION_REF_SEQUENCE_DICTIONARY = _dir_seq_extraction / 'fasta_ref.io'
OUTPUT_SEQ_EXTRACTION_HAPLOTYPECALLER = _dir_seq_extraction / 'haplotypecaller' / 'vcf.io'
OUTPUT_SEQ_EXTRACTION_VARIANTFILTRATION = _dir_seq_extraction / 'variantfiltration' / 'vcf.io'
OUTPUT_SEQ_EXTRACTION_SELECTVARIANTS = _dir_seq_extraction / 'selectvariants' / 'vcf.io'
OUTPUT_SEQ_EXTRACTION_SAMTOOLS_DEPTH = _dir_seq_extraction / 'samtools_depth' / 'tsv.io'
OUTPUT_SEQ_EXTRACTION_SAMTOOLS_DEPTH_STATS = _dir_seq_extraction / 'samtools_depth_stats' / 'samtools_depth_informs.io'
OUTPUT_SEQ_EXTRACTION_VCF_INDEL_SCAN = _dir_seq_extraction / 'vcf_indel_scan' / 'vcf_indel_scan_informs.io'
OUTPUT_SEQ_EXTRACTION_REGION_CALCULATOR_INFORMS = _dir_seq_extraction / 'region_calculator' / 'region_calculator_informs.io'
OUTPUT_SEQ_EXTRACTION_REGION_CALCULATOR_INTERVALS = _dir_seq_extraction / 'region_calculator' / 'tsv.io'
OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE = _dir_seq_extraction / 'consensus_sequence' / 'fasta.io'
OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE_INDEX_PREFIX = _dir_seq_extraction / 'consensus_sequence' / 'index_genome_prefix.io'
OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE_ITERATIVE = _dir_seq_extraction / 'seq_extraction_iterative' / 'fasta.io'
