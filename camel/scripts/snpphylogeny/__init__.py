import os

_current_dir = os.path.dirname(__file__)
SNAKEFILE_SAMTOOLS_CALLING_ALL = os.path.join(_current_dir, 'snakefile', 'samtools_calling_all.smk')
SNAKEFILE_SAMTOOLS_FILTERING_ALL = os.path.join(_current_dir, 'snakefile', 'samtools_filtering_all.smk')
