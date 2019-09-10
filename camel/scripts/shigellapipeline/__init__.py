from pathlib import Path

_current_dir = Path(__file__).parent
CITATIONS_HTML = str(_current_dir / 'citations.html')
CONFIG_TEMPLATE = str(_current_dir / 'config' / 'config_template.yml')
CONFIG_DATA = str(_current_dir / 'config' / 'config_data.yml')
SNAKEFILE_MAIN = str(_current_dir / 'snakefile' / 'main.smk')
SNAKEFILE_SUBSPECIES = str(_current_dir / 'snakefile' / 'subspecies_identification.smk')
SNAKEFILE_FLEXNERI_TYPE = str(_current_dir / 'snakefile' / 'flexneritype.smk')
