import os

_current_dir = os.path.dirname(__file__)
SNAKEFILE_MAIN = os.path.join(_current_dir, 'snakefile', 'main.smk')
CONFIG_TEMPLATE = os.path.join(_current_dir, 'config', 'config_template.yml')
CONFIG_DATA = os.path.join(_current_dir, 'config', 'config_data.yml')
CITATIONS_HTML = os.path.join(_current_dir, 'static', 'citations.html')
