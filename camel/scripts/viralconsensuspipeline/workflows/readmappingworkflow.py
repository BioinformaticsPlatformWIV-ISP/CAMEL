import dataclasses
import json
from pathlib import Path
from typing import Optional, Any, Dict, List

import pkg_resources

from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.loggers import logger
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils


@dataclasses.dataclass
class ReadMappingOutput:
    """
    Holder for the output of the read mapping workflow.
    """
    path_bam: Path
    path_bed_low_cov: Path
    stats: Dict[str, Any]
    informs: List[Dict[str, Any]]


class ReadMappingWorkflow(object):
    """
    Maps reads to an input FASTA file and extracts depth statistics.
    Illumina data: BWA (separately mapping SE and PE reads and merges them afterward).
    ONT data: Minimap2
    """

    def __init__(self, dir_: Path) -> None:
        """
        Initializes this workflow.
        :param dir_: Working directory
        :return: None
        """
        self._dir = dir_
        if not self._dir.exists():
            logger.info(f'Creating working directory: {self._dir}')
            self._dir.mkdir(parents=True)

    def run(self, fastq_in: FastqInput, fasta_ref: Path, threads: int = 8, prefix: Optional[str] = 'mapping',
            gap_len_cutoff: int = 6, gap_depth_cutoff: int = 10) -> ReadMappingOutput:
        """
        Runs the read mapping workflow.
        :param fastq_in: Input FASTQ data
        :param fasta_ref: Reference FASTA file
        :param prefix: Prefix for the output files
        :param threads: Number of threads
        :param gap_len_cutoff: Min gap length
        :param gap_depth_cutoff: Min gap coverage
        :return: Output holder
        """
        # Create config file
        path_config = SnakePipelineUtils.generate_config_file({
            'input': {
                'prefix': prefix,
                'FASTA': str(fasta_ref)
            },
            'mapper': 'bwa' if fastq_in.read_type == 'illumina' else 'minimap2',
            'low_depth': {'gap_depth_cutoff': gap_depth_cutoff, 'gap_len_cutoff': gap_len_cutoff}
        }, self._dir)

        # Create input directory
        dir_input = self._dir / 'input'
        dir_input.mkdir(exist_ok=True)
        SnakemakeUtils.dump_object(fastq_in.to_fq_dict(), dir_input / 'fq_dict.io')

        # Run Snakefile
        path_snakefile = pkg_resources.resource_filename(
            'camel', 'scripts/viralconsensuspipeline/workflows/readmappingworkflow.smk')
        SnakePipelineUtils.run_snakemake(path_snakefile, str(path_config), [], working_dir=self._dir, threads=threads)

        # Collect output
        path_bam = SnakemakeUtils.load_object(self._dir / 'sam_to_bam' / 'bam.io')[0].path
        path_bed = SnakemakeUtils.load_object(self._dir / 'low_depth' / 'bed.io')[0].path
        with (self._dir / 'stats_mapping.json').open() as handle:
            stats = json.load(handle)
        with (self._dir / 'informs_all.json').open() as handle:
            informs = json.load(handle)
        return ReadMappingOutput(path_bam=path_bam, path_bed_low_cov=path_bed, stats=stats, informs=informs)
