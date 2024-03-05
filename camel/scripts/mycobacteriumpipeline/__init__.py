import os

_current_dir = os.path.dirname(__file__)
SNAKEFILE_51SNP = os.path.join(_current_dir, 'snakefile', 'assay51snp.smk')
SNAKEFILE_AMR = os.path.join(_current_dir, 'snakefile', 'amr.smk')
SNAKEFILE_MAIN = os.path.join(_current_dir, 'snakefile', 'main.smk')
SNAKEFILE_HSP65 = os.path.join(_current_dir, 'snakefile', 'hsp65.smk')
SNAKEFILE_CSB_RD = os.path.join(_current_dir, 'snakefile', 'csb_rd.smk')
SNAKEFILE_SNP_LINEAGE = os.path.join(_current_dir, 'snakefile', 'snplineage.smk')
SNAKEFILE_SNPIT = os.path.join(_current_dir, 'snakefile', 'snpit.smk')
SNAKEFILE_SPOLIGOTYPING = os.path.join(_current_dir, 'snakefile', 'spoligotyping.smk')
CONFIG_TEMPLATE = os.path.join(_current_dir, 'config', 'config_template.yml')
CONFIG_DATA = os.path.join(_current_dir, 'config', 'config_data.yml')
AMR_CIRCOS_TEMPLATE = os.path.join(_current_dir, 'static', 'amr_circos_template.txt')
