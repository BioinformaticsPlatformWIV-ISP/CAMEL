from pathlib import Path

_current_dir = Path(__file__).parent
CONFIG_DATA = str(_current_dir / 'config' / 'config_data.yml')
SNAKEFILE_MAIN = str(_current_dir / 'snakefile' / 'main.smk')
SNAKEFILE_SHIGATYPER = str(_current_dir / 'snakefile' / 'shigatyper.smk')
SNAKEFILE_SHIGEIFINDER = str(_current_dir / 'snakefile' / 'shigeifinder.smk')
