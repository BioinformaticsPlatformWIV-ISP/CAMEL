import os
_config_folder = os.path.dirname(os.path.realpath(__file__))

DB_CONFIG = os.path.join(_config_folder, 'db.yml')
LOGGING_CONFIG = os.path.join(_config_folder, 'logging.yml')
