import os

__CURRENT_DIR = os.path.dirname(__file__)

WORKFLOW_ASSEMBLY = os.path.join(__CURRENT_DIR, 'assembly.snakefile')
WORKFLOW_GENE_DETECTION = os.path.join(__CURRENT_DIR, 'gene_detection.snakefile')
WORKFLOW_READ_TRIMMING = os.path.join(__CURRENT_DIR, 'read_trimming.snakefile')
