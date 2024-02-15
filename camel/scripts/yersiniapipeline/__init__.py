import os

_current_dir = os.path.dirname(__file__)
CONFIG_DATA = os.path.join(_current_dir, 'config', 'config_data.yml')
SNAKEFILE_MAIN = os.path.join(_current_dir, 'snakefile', 'main.smk')
#TODO: serogroup? other custom snakefiles?