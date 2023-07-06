from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Any, Optional

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import short_read_polishing


class ShortReadPolishingWrapper(object):
    """
    This class is used as a wrapper class around the Illumina short read polishing Snakemake workflow.
    """

    @dataclass
    class PolishingOutput:
        fasta_contigs: Path
        log_file: Optional[Path] = None

    def __init__(self, working_dir: Path) -> None:
        """
        Initializes the read trimming wrapper.
        :param working_dir: Working directory
        """
        self._working_dir = Path(working_dir)
        self._output = None

    def run_workflow(self, pe_reads: Path, reference: Path = None, threads: int = 8) -> None:
        """
        Runs the read polishing workflow.
        :param pe_reads: Input PE FASTQ reads
        :param reference: genome to polish
        :param threads: Number of threads to use
        :return: None
        """
        # Create config file
        config_data = {
            'working_dir': str(self._working_dir),
            'polishing': {'polypolish': [], 'polca': []},
        }
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._working_dir)

        # Dump the input files in an IO file
        io_pickle_in = Path(self._working_dir / short_read_polishing.INPUT_READS_FASTQ)
        io_pickle_in.parent.mkdir(exist_ok=True, parents=True)
        SnakemakeUtils.dump_object([ToolIOFile(pe_reads)], io_pickle_in)

        # Dump the reference file in an IO file
        io_pickle_fasta_in = Path(self._working_dir / short_read_polishing.INPUT_ASSEMBLY_FASTA)
        io_pickle_fasta_in.parent.mkdir(exist_ok=True, parents=True)
        SnakemakeUtils.dump_object([ToolIOFile(reference)], io_pickle_in)

        # Output files
        output_files = {
            'FASTA': self._working_dir / short_read_polishing.OUTPUT_POLISHING_FASTA
        }
        SnakePipelineUtils.run_snakemake(
            short_read_polishing.SNAKEFILE_POLISHING, config_file, list(output_files.values()), self._working_dir,
            threads)
        self.__set_output(output_files)

    def __set_output(self, output_files: Dict[str, Path]) -> None:
        """
        Sets the output of this tool.
        :param output_files: Output files by key.
        :return: None
        """
        log_path = self._working_dir / 'camel.log'
        self._output = ShortReadPolishingWrapper.PolishingOutput(
            fasta_contigs=output_files['FASTA'],
            log_file=log_path if log_path.exists() else None
        )

    @property
    def output(self) -> 'ShortReadPolishingWrapper.PolishingOutput':
        """
        Returns the report section generated during the read trimming pipeline.
        :return: Trimming section
        """
        return self._output
