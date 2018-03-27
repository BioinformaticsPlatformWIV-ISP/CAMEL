import os

_current_dir = os.path.abspath(os.path.dirname(__file__))

SNAKEFILE_VTEC_MAIN = os.path.join(_current_dir, 'main.snakefile')
