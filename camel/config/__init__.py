from pathlib import Path

_config_folder = Path(__file__).parent

LOGGING_CONFIG = _config_folder / 'logging.yml'
MAIN_CONFIG = _config_folder / 'main.yml'
