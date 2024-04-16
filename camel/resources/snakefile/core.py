from pathlib import Path
from typing import Dict, Any

from camel.resources.snakefile import trimming_illumina, trimming_ont, human_read_scrubbing

SNAKEFILE_CORE = f'{Path(__file__).parent / Path(__file__).stem}.smk'
INPUT_FASTA_IO = Path('input', 'fasta.io')
OUTPUT_TSV_SUMMARY_INIT = Path('summary', 'summary-init.tsv')
OUTPUT_HTML_CITATIONS = Path('report', 'html-citations.io')


def get_fastq_input_downsampling(config: Dict[str, Any], read_key: str) -> Path:
    """
    Returns the fastq input for the downsampling step.
    :param config: Snakemake configuration
    :param read_key: fastq_se or fastq_pe
    :return: Path to the input file
    """
    if 'human_read_scrubbing' in config.get('analyses', []):
        return Path(config['working_dir']) / str(human_read_scrubbing.OUTPUT_SCRUBBING_FASTQ).format(input_format=read_key)
    else:
        return Path(config['working_dir']) / str(human_read_scrubbing.INPUT_SCRUBBING_FASTQ).format(input_format=read_key)


def get_fq_input(input_type: str, dir_: Path):
    """
    Returns the FQ input (which is unpacked into the Snakemake input).
    :param input_type: Input type
    :param dir_: Working directory
    :return: Dictionary with strings as keys and paths as values
    """
    if input_type in ('fasta', 'fasta_vcf'):
        return {}
    elif input_type == 'illumina':
        return {'FASTQ_PE': dir_ / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_DICT}
    elif input_type == 'ont':
        return {'FASTQ_SE': dir_ / trimming_ont.OUTPUT_TRIMMING_ONT_DICT}
    elif input_type == 'hybrid':
        return {
            'FASTQ_PE': dir_ / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_DICT,
            'FASTQ_SE': dir_ / trimming_ont.OUTPUT_TRIMMING_ONT_DICT
        }
    else:
        raise ValueError(f'Invalid input type: {input_type}')
