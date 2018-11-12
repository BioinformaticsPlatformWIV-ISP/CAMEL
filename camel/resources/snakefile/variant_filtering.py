from typing import Dict, Any

import os

OUTPUT_VARIANT_FILTERING_STATS = os.path.join('variant_filtering', 'stats', 'json.io')
OUTPUT_VARIANT_FILTERING_SUMMARY = os.path.join('variant_filtering', 'summary.tsv')
OUTPUT_VARIANT_FILTERING_VCF = os.path.join('variant_filtering', 'unzip', 'vcf.io')
OUTPUT_VARIANT_FILTERING_VCF_GZ = os.path.join('variant_filtering', 'zscore', 'vcf_gz.io')


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
