import os

_current_dir = os.path.dirname(__file__)
CITATIONS_HTML = os.path.join(_current_dir, 'citations.html')
CONFIG_TEMPLATE = os.path.join(_current_dir, 'config', 'config_template.yml')
CONFIG_DATA = os.path.join(_current_dir, 'config', 'config_data.yml')
SNAKEFILE_MAIN = os.path.join(_current_dir, 'snakefile', 'main.smk')
SNAKEFILE_SEROTYPE = os.path.join(_current_dir, 'snakefile', 'serotypedetection.smk')
