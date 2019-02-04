import os

_current_dir = os.path.dirname(__file__)
SNAKEFILE_MAIN = os.path.join(_current_dir, 'snakefile', 'main.smk')
SNAKEFILE_SEROTYPE = os.path.join(_current_dir, 'snakefile', 'serotype_detection.smk')
CONFIG_TEMPLATE = os.path.join(_current_dir, 'config', 'config_template.yml')
CITATIONS_HTML = os.path.join(_current_dir, 'citations.html')
