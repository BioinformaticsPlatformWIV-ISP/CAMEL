from pathlib import Path

SNAKEFILE_VARIANT_CALLING = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_variant_calling = Path('variant_calling')
INPUT_VARIANT_CALLING_BAM = _dir_variant_calling / 'read_mapping' / 'bam.io'
OUTPUT_VARIANT_CALLING_BAM = _dir_variant_calling / 'alignment_sorting' / 'bam-sorted.io'
OUTPUT_VARIANT_CALLING_BAM_INDEX = _dir_variant_calling / 'alignment_sorting' / 'bam-index.io'
OUTPUT_VARIANT_CALLING_FASTA = _dir_variant_calling / 'reference' / 'fasta.io'
OUTPUT_VARIANT_CALLING_FASTA_INDEX = _dir_variant_calling / 'reference' / 'fasta-index.io'
OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ = _dir_variant_calling / 'norm' / 'vcf_gz-indexed.io'
OUTPUT_VARIANT_CALLING_UNFILTERED_VCF = _dir_variant_calling / 'unzip_vcf' / 'vcf.io'