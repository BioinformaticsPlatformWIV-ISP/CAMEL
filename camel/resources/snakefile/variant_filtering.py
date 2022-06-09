from pathlib import Path
from typing import Dict, Any

SNAKEFILE_VARIANT_FILTERING = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_variant_filtering = Path('variant_filtering')

OUTPUT_VARIANT_FILTERING_STATS = _dir_variant_filtering / 'stats' / 'json.io'
OUTPUT_VARIANT_FILTERING_INFORMS_ALL = _dir_variant_filtering / 'informs_all.io'
OUTPUT_VARIANT_FILTERING_SUMMARY = _dir_variant_filtering / 'summary.tsv'
OUTPUT_VARIANT_FILTERING_VCF = _dir_variant_filtering / 'regions' / 'vcf.io'
OUTPUT_VARIANT_FILTERING_VCF_GZ = _dir_variant_filtering / 'zscore' / 'vcf_gz.io'


def get_filtering_param(config: Dict[str, Any], filter_key: str, param_name: str) -> Any:
    """
    Returns the value of the given filtering parameter, returns 'None' when the parameter is not set.
    :param config: Snakemake configuration
    :param filter_key: Filter key
    :param param_name: Parameter name
    :return: Parameter value if set, None otherwise
    """
    if filter_key not in config['variant_filtering']:
        return None
    return config['variant_filtering'][filter_key].get(param_name)


def filter_is_disabled(filter_key: str, config: Dict[str, Any]) -> bool:
    """
    Returns True if the given filter is disabled, False otherwise.
    :param filter_key: Filter key
    :param config: Snakemake config
    :return: True if disabled, false otherwise
    """
    if filter_key not in config['variant_filtering']:
        return False
    return config['variant_filtering'][filter_key].get('disabled', False)
