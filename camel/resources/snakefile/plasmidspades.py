from pathlib import Path

import snakemake

from camel.app.snakemake.snakemakeutils import SnakemakeUtils

SNAKEFILE_PLASMID_SPADES = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_plasmidspades = Path('plasmidspades')

OUTPUT_PLASMIDSPADES_REPORT = _dir_plasmidspades / 'report' / 'html.io'
OUTPUT_PLASMIDSPADES_REPORT_EMPTY = _dir_plasmidspades / 'report' / 'html_empty.io'
OUTPUT_PLASMIDSPADES_GENE_DETECTION_REPORT = _dir_plasmidspades / 'gene_detection' / '{db}' / 'html.io'
OUTPUT_PLASMIDSPADES_GENE_DETECTION_REPORT_EMPTY = _dir_plasmidspades / 'gene_detection' / '{db}' / 'html-empty.io'
OUTPUT_PLASMIDSPADES_INFORMS = _dir_plasmidspades / 'informs.io'
OUTPUT_PLASMIDSPADES_SUMMARY = _dir_plasmidspades / 'summary.tsv'


def plasmidspades_successful(checkpoint: snakemake.workflow.Checkpoint) -> bool:
    """
    Function to check if the plasmidSPAdes assembly was successful.
    :param checkpoint: Checkpoint
    :return: True if successful, False otherwise
    """
    return len(SnakemakeUtils.load_object(Path(checkpoint.get().output.FASTA_Contig))) > 0
