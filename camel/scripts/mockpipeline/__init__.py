from importlib.resources import files
from pathlib import Path

SNAKEFILE_MAIN = Path(str(files('camel').joinpath('scripts/mockpipeline/main.smk')))
CONFIG_DATA = Path(str(files('camel').joinpath('scripts/mockpipeline/config_data.yml')))
