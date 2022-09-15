from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional, Dict, Any

import operator

from camel.resources.snakefile.assembly_spades import OUTPUT_ASSEMBLY_MAPPING_INFORMS, OUTPUT_ASSEMBLY_DEPTH_INFORMS
from camel.resources.snakefile.variant_calling import OUTPUT_VARIANT_CALLING_MAPPING_INFORMS, \
    OUTPUT_VARIANT_CALLING_DEPTH_INFORMS

SNAKEFILE_QUALITY_CHECKS = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_QUALITY_CHECKS_REPORT = Path('quality_checks') / 'report' / 'html.io'
OUTPUT_QUALITY_CHECKS_SUMMARY = Path('quality_checks') / 'summary' / 'summary_out.tsv'
OUTPUT_QUALITY_CHECKS_REPORT_JSON = Path('quality_checks') / 'report' / 'informs.json'


def get_mapping_rate_informs(config) -> Path:
    """
    Returns the input for the mapping rate.
    :param config: Snakemake config
    :return: Mapping rate informs
    """
    mode = config['quality_checks'].get('coverage_mode', 'assembly')
    if mode == 'assembly':
        return Path(config['working_dir']) / OUTPUT_ASSEMBLY_MAPPING_INFORMS
    elif mode == 'ref':
        return Path(config['working_dir']) / OUTPUT_VARIANT_CALLING_MAPPING_INFORMS
    raise ValueError(f"Invalid coverage mode: '{mode}'")


def get_depth_informs(config) -> Path:
    """
    Returns the input for the depth.
    :param config: Snakemake config
    :return: Depth informs
    """
    mode = config['quality_checks'].get('coverage_mode', 'assembly')
    if mode == 'assembly':
        return Path(config['working_dir']) / OUTPUT_ASSEMBLY_DEPTH_INFORMS
    elif mode == 'ref':
        return Path(config['working_dir']) / OUTPUT_VARIANT_CALLING_DEPTH_INFORMS
    raise ValueError(f"Invalid coverage mode: '{mode}'")


@dataclass
class QCCheck:
    key: str
    full_name: str
    threshold_warn: float
    threshold_fail: float
    fmt_string_value: Optional[str]
    value_should_exceed: Optional[bool] = True

    def to_dict(self, value: Optional[float] = None) -> Dict[Any, str]:
        """
        Converts the QC check to a dictionary.
        :param value: Metric value
        :return: Dictionary
        """
        data = asdict(self)
        comp = operator.lt if self.value_should_exceed else operator.gt
        if value is None:
            data['status'] = 'Skipped'
        elif comp(value, self.threshold_fail):
            data['status'] = 'Failed'
        elif comp(value, self.threshold_warn):
            data['status'] = 'Warning'
        else:
            data['status'] = 'OK'
        data['value'] = value
        return data


QC_CHECKS_BY_KEY = {qc.key: qc for qc in [
    QCCheck('cgmlst', 'Typing loci detected (%)', 95.0, 90.0, '{:.2f}%'),
    QCCheck('kraken', 'Kraken: contaminants', 1.0, 5.0, '{:.2f}%', False),
    QCCheck('map_rate_ref', 'Reads mapping to reference genome', 95.0, 90.0, '{:.2f}%'),
    QCCheck('cov_ref', 'Coverage against reference genome', 20, 10, '{:.2f}x'),
    QCCheck('map_rate_assembly', 'Reads mapping to the assembled contigs', 95.0, 90.0, '{:.2f}%'),
    QCCheck('cov_assembly', 'Coverage against the assembled contigs', 20, 10, '{:.2f}x'),
    QCCheck('fqc_gc', 'FastQC: GC-content deviation', 2.0, 4.0, '{:.2f}%', False),
    QCCheck('fqc_avg_qual', 'FastQC: Average quality score', 30.0, 25.0, '{:.0f}'),
    QCCheck('fqc_n_fraction', 'FastQC: Max. N-fraction', 0.005, 0.010, '{:.4f}', False),
    QCCheck('fqc_per_base', 'FastQC: Per-base sequence content', 3.0, 6.0, '{:.2f}%', False),
    QCCheck('fqc_qscore', 'FastQC: Q-score drop', 200, 150, '{:.0f}'),
    QCCheck('fqc_seq_len', 'FastQC: Sequence length distribution', 66.67, 40.00, '{:.2f}%')
]}
