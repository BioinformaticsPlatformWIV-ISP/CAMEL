import os

_current_dir = os.path.dirname(__file__)
SNAKEFILE_MAIN = os.path.join(_current_dir, 'snakefile', 'main.smk')
SNAKEFILE_SEROGROUP_DETERMINATION = os.path.join(_current_dir, 'snakefile', 'serogroup_determination.smk')
CONFIG_TEMPLATE = os.path.join(_current_dir, 'config', 'config_template.yml')
CITATIONS_HTML = os.path.join(_current_dir, 'citations.html')
