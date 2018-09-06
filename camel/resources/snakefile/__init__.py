import os
_current_dir = os.path.dirname(os.path.realpath(__file__))

SNAKEFILE_GENE_DETECTION = os.path.join(_current_dir, 'gene_detection.smk')
SNAKEFILE_SEQUENCE_TYPING = os.path.join(_current_dir, 'sequence_typing.smk')
