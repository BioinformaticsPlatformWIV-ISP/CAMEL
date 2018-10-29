import os
_current_dir = os.path.dirname(os.path.realpath(__file__))

SNAKEFILE_ASSEMBLY_SPADES = os.path.join(_current_dir, 'assembly_spades.smk')
SNAKEFILE_GENE_DETECTION = os.path.join(_current_dir, 'gene_detection.smk')
SNAKEFILE_READ_TRIMMING = os.path.join(_current_dir, 'read_trimming.smk')
SNAKEFILE_SEQUENCE_TYPING = os.path.join(_current_dir, 'sequence_typing.smk')
SNAKEFILE_VARIANT_CALLING = os.path.join(_current_dir, 'variant_calling.smk')
SNAKEFILE_VARIANT_FILTERING = os.path.join(_current_dir, 'variant_filtering.smk')
