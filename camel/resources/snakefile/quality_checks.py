import operator
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional, Dict, Any

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
    explanation: Optional[str] = None

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
    # General
    QCCheck('cgmlst', 'Typing loci detected (%)', 95.0, 90.0, '{:.2f}%'),
    QCCheck('kraken', 'Kraken: contaminants', 1.0, 5.0, '{:.2f}%', False),
    QCCheck('map_rate_ref', 'Reads mapping to reference genome', 95.0, 90.0, '{:.2f}%'),
    QCCheck('cov_ref', 'Coverage against reference genome', 20, 10, '{:.2f}x'),
    QCCheck('map_rate_assembly', 'Reads mapping to the assembled contigs', 95.0, 90.0, '{:.2f}%'),
    QCCheck('cov_assembly', 'Coverage against the assembled contigs', 20, 10, '{:.2f}x'),

    # Illumina
    QCCheck('fqc_gc', 'FastQC: GC-content deviation', 2.0, 4.0, '{:.2f}%', False,
            explanation='checks if the detected GC content is close enough to the expected GC content for this organism'
                        ' (<b>{:.2f}%</b>).'),
    QCCheck('fqc_avg_qual', 'FastQC: Average quality score', 30.0, 25.0, '{:.0f}',
            explanation='checks if the average read quality is above the given threshold.'),
    QCCheck('fqc_n_fraction', 'FastQC: Max. N-fraction', 0.005, 0.010, '{:.4f}', False,
            explanation='checks if the maximal N fraction at any read position is below the given threshold.'),
    QCCheck('fqc_per_base', 'FastQC: Per-base sequence content', 3.0, 6.0, '{:.2f}%', False,
            explanation='checks if the difference between A-T and C-G is below the given threshold at every position. '
                        'The first 20 and last 5 bases of the reads are skipped, as the peaks there can be caused by '
                        'the library kit or trimming artifacts.'),
    QCCheck('fqc_qscore', 'FastQC: Q-score drop', 200, 150, '{:.0f}',
            explanation='checks whether the average position in the reads where the mean Q-score drops below <b>30</b> '
                        'is above the given threshold.'),
    QCCheck('fqc_seq_len', 'FastQC: Sequence length distribution', 66.67, 40.00, '{:.2f}%',
            explanation='checks if the median read length of the trimmed reads is below a threshold compared to the '
                        'mode length of the raw input reads (<b>{}</b>).'),

    # Nanopore
    QCCheck('nanoplot_len', 'NanoPlot: Median read length', 500, 250, '{:.0f}'),
    QCCheck('nanoplot_qual', 'NanoPlot: Median read quality', 10, 8, '{:.2f}%'),
    QCCheck('seqkit_gc', 'seqkit: GC-content deviation', 10, 8, '{:.2f}%',
            explanation='checks if the detected GC content is close enough to the expected GC content for this organism'
                        ' (<b>{:.2f}%</b>).')
]}
