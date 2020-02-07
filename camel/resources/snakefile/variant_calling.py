from pathlib import Path


SNAKEFILE_VARIANT_CALLING = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_variant_calling = Path('variant_calling')

OUTPUT_VARIANT_CALLING_REPORT = _dir_variant_calling / 'report' / 'html.io'
OUTPUT_VARIANT_CALLING_SUMMARY = _dir_variant_calling / 'summary' / 'summary_out.tsv'
OUTPUT_VARIANT_CALLING_BAM = _dir_variant_calling / 'alignment_sorting' / 'bam-sorted.io'
OUTPUT_VARIANT_CALLING_DEPTH_INFORMS = _dir_variant_calling / 'depth' / 'informs.io'
OUTPUT_VARIANT_CALLING_DEPTH_TSV = _dir_variant_calling / 'depth' / 'tsv.io'
OUTPUT_VARIANT_CALLING_UNFILTERED_VCF = _dir_variant_calling / 'unzip_vcf' / 'vcf.io'
OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ = _dir_variant_calling / 'calling' / 'vcf_gz_indexed.io'
OUTPUT_VARIANT_CALLING_CONSENSUS = _dir_variant_calling / 'consensus' / 'fasta.io'
OUTPUT_VARIANT_CALLING_FILTERED_VCF = _dir_variant_calling / 'unzip_vcf_filtered' / 'vcf.io'
OUTPUT_VARIANT_CALLING_FILTERED_VCF_GZ = _dir_variant_calling / 'filter_zscore' / 'vcf_gz-indexed.io'
OUTPUT_VARIANT_CALLING_MAPPING_INFORMS = _dir_variant_calling / 'read_mapping' / 'informs.io'
OUTPUT_VARIANT_CALLING_INFORMS_ALL = _dir_variant_calling / 'informs_all.io'
