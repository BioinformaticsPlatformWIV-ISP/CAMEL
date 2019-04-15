import os

OUTPUT_VARIANT_CALLING_REPORT = os.path.join('variant_calling', 'report', 'html.io')
OUTPUT_VARIANT_CALLING_SUMMARY = os.path.join('variant_calling', 'summary', 'summary_out.tsv')
OUTPUT_VARIANT_CALLING_BAM = os.path.join('variant_calling', 'alignment_sorting', 'bam-sorted.io')
OUTPUT_VARIANT_CALLING_UNFILTERED_VCF = os.path.join('variant_calling', 'unzip_vcf', 'vcf.io')
OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ = os.path.join('variant_calling', 'calling', 'vcf_gz_indexed.io')
OUTPUT_VARIANT_CALLING_CONSENSUS = os.path.join('variant_calling', 'consensus', 'fasta.io')
OUTPUT_VARIANT_CALLING_FILTERED_VCF = os.path.join('variant_calling', 'unzip_vcf_filtered', 'vcf.io')
OUTPUT_VARIANT_CALLING_FILTERED_VCF_GZ = os.path.join('variant_calling', 'filter_zscore', 'vcf_gz-indexed.io')
OUTPUT_VARIANT_CALLING_MAPPING_INFORMS = os.path.join('variant_calling', 'read_mapping', 'informs.io')
OUTPUT_VARIANT_CALLING_INFORMS_ALL = os.path.join('variant_calling', 'commands.txt')
