from pathlib import Path

FOLDER_TRIMMING_ILLUMINA = Path('trimming_illumina')
SNAKEFILE_TRIMMING_ILLUMINA = f'{Path(__file__).parent / Path(__file__).stem}.smk'

INPUT_TRIMMING_FASTQ = FOLDER_TRIMMING_ILLUMINA / 'input' / 'fastq.io'

# Report and summary
OUTPUT_TRIMMING_ILLUMINA_REPORT = FOLDER_TRIMMING_ILLUMINA / 'report' / 'html.io'
OUTPUT_TRIMMING_ILLUMINA_SUMMARY = FOLDER_TRIMMING_ILLUMINA / 'summary' / 'summary_trim.tsv'

# Trimming
OUTPUT_TRIMMING_ILLUMINA_READS_PE = FOLDER_TRIMMING_ILLUMINA / 'reads' / 'fastq-pe.io'
OUTPUT_TRIMMING_ILLUMINA_READS_SE_FWD = FOLDER_TRIMMING_ILLUMINA / 'reads' / 'fastq-se-fwd.io'
OUTPUT_TRIMMING_ILLUMINA_READS_SE_REV = FOLDER_TRIMMING_ILLUMINA / 'reads' / 'fastq-se-rev.io'
OUTPUT_TRIMMING_ILLUMINA_DICT = FOLDER_TRIMMING_ILLUMINA / 'reads' / 'fq_dict.io'
OUTPUT_TRIMMING_ILLUMINA_INFORMS = FOLDER_TRIMMING_ILLUMINA / 'report' / 'informs.io'

# FastQC
OUTPUT_TRIMMING_ILLUMINA_FASTQC_TXT_PRE = FOLDER_TRIMMING_ILLUMINA / 'fastqc-pre' / 'txt.io'
OUTPUT_TRIMMING_ILLUMINA_FASTQC_HTML_PRE = FOLDER_TRIMMING_ILLUMINA / 'fastqc-pre' / 'html.io'
OUTPUT_TRIMMING_ILLUMINA_FASTQC_TXT_POST = FOLDER_TRIMMING_ILLUMINA / 'fastqc-post' / 'txt.io'
OUTPUT_TRIMMING_ILLUMINA_FASTQC_HTML_POST = FOLDER_TRIMMING_ILLUMINA / 'fastqc-post' / 'html.io'


def select_fastq_output(config: dict) -> Path:
    """
    Selects the PE FASTQ output based on the Snakemake configuration.
    :param config: Config data
    :return: Path to PE FASTQ output
    """
    method = config['read_trimming'].get('method', 'trimmomatic')
    if method == 'trimmomatic':
        return Path('trimming_illumina') / 'trimmomatic' / 'fastq-pe.io'
    else:
        return Path('trimming_illumina') / 'fastp' / 'fastq-pe.io'

def select_informs(config: dict) -> Path:
    """
    Selects the trimming informs based on the Snakemake configuration.
    :param config: Config data
    :return: Path to the trimming informs
    """
    method = config['read_trimming'].get('method', 'trimmomatic')
    if method == 'trimmomatic':
        return Path('trimming_illumina') / 'trimmomatic' / 'informs.io'
    else:
        return Path('trimming_illumina') / 'fastp' / 'informs.io'

def select_report_output(config: dict) -> Path:
    """
    Selects the trimming report based on the Snakemake configuration.
    :param config: Config data
    :return: Path to the output report
    """
    method = config['read_trimming'].get('method', 'trimmomatic')
    if method == 'trimmomatic':
        return Path('trimming_illumina') / 'report' / 'trimmomatic' / 'html.io'
    else:
        return Path('trimming_illumina') / 'report' / 'fastp' / 'html.io'
