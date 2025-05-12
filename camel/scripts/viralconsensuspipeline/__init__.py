from pathlib import Path

import pkg_resources

SNAKEFILE_MAIN = Path(pkg_resources.resource_filename('camel', 'scripts/viralconsensuspipeline/snakefile/main.smk'))
