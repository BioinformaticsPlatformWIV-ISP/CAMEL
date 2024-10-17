from pathlib import Path
from collections import defaultdict
from typing import List, Tuple, Union

_current_dir = Path(__file__).parent
CONFIG_DATA = _current_dir / 'config' / 'config_data.yml'
SNAKEFILE_MAIN = _current_dir / 'snakefile' / 'main.smk'


def add_content_spifinder(
        structure: List[Tuple], input_type: str, reports_spifinder: List[Union[Path, str]]) -> None:
    """
    Adds the report content for SPIFinder.
    :param structure: Report structure
    :param input_type: Input type
    :param reports_spifinder: SPIFinder output reports
    :return: None
    """
    # Create dictionaries with the technology as key and the reports as values
    report_ds_by_input_format = {
        p_html.parent.name: p_html for p_html in [Path(x) for x in reports_spifinder]}

    # Add the report content
    if input_type == 'fasta':
        structure.append(('Pathogenicity island determination', 'spif', [report_ds_by_input_format['fasta']]))
    elif input_type == 'illumina':
        structure.append(('Pathogenicity island determination', 'spif', [report_ds_by_input_format['fastq_pe']]))
    elif input_type == 'ont':
        structure.append(('Pathogenicity island determination', 'spif', [report_ds_by_input_format['fastq_se']]))
    elif input_type == 'hybrid':
        structure.append(
            ('Pathogenicity island determination - Illumina', 'spif_ilmn', [report_ds_by_input_format['fastq_pe']]))
        structure.append(
            ('Pathogenicity island determination - ONT', 'spif_ont', [report_ds_by_input_format['fastq_se']]))


def add_content_serotyping_salmonella(
        structure: List[Tuple], input_type: str, reports_serotyping: List[Union[Path, str]]) -> None:
    """
    Adds the report content for the Salmonella serotyping assays SISTR and SeqSero2.
    :param structure: Report structure
    :param input_type: Input type
    :param reports_serotyping: Salmonella serotyping output reports
    :return: None
    """
    report_ds_by_input_format = defaultdict(list)

    for report_path in reports_serotyping:
        p_html = Path(report_path)
        parent_name = p_html.parent.name
        report_ds_by_input_format[parent_name].append(p_html)

    # Add the report content
    if input_type == 'fasta':
        structure.append(('Serotyping', 'sero', [Path(x) for x in report_ds_by_input_format['fasta']]))
    elif input_type == 'illumina':
        structure.append(('Serotyping', 'sero', [Path(x) for x in report_ds_by_input_format['fastq_pe']]))
    elif input_type == 'ont':
        structure.append(('Serotyping', 'sero', [Path(x) for x in report_ds_by_input_format['fastq_se']]))
    elif input_type == 'hybrid':
        structure.append(
            ('Serotyping - Illumina', 'sero_ilmn', [Path(x) for x in report_ds_by_input_format['fastq_pe']]))
        structure.append(('Serotyping - ONT', 'sero_ont', [Path(x) for x in report_ds_by_input_format['fastq_se']]))
