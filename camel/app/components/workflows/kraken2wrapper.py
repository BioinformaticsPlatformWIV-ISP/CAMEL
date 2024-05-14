from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Dict, Any, Union

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import contamination_check_kraken


@dataclass
class Kraken2Output:
    report_section: HtmlReportSection
    tsv_summary: Path
    informs: Dict[str, Any]
    informs_commands: List[Dict[str, Any]]
    log_file: Optional[Path]


class Kraken2Wrapper(object):
    """
    This class is used as a wrapper class around the Kraken2 contamination check workflow.
    """

    def __init__(self, working_dir: Union[str, Path]) -> None:
        """
        Initializes the read trimming helper.
        :param working_dir: Working directory
        :return: None
        """
        self._working_dir = Path(working_dir)
        self._output = None

    def run_workflow(
            self, sample_name: str, input_path: Union[Path, List[Path]], input_type: str, expected_species: str, db: Path, level_of_depth: str = 'S',
            threads: int = 8) -> None:
        """
        Runs the kraken2 workflow.
        :param input_path: Path to the input
        :param input_type: input type
        :param sample_name: Sample name
        :param expected_species: Expected species
        :param db: Database
        :param level_of_depth: Species ('S') or Genus ('G') level of contamination check
        :param threads: Number of threads
        :return: None
        """
        if not self._working_dir.exists():
            self._working_dir.mkdir(parents=True)
        if input_type == 'illumina':
            fastq_input = FastqInput('illumina', pe=[ToolIOFile(x) for x in input_path])
            SnakemakeUtils.dump_object(fastq_input.to_fq_dict(), self._working_dir / 'fq_dict.io')
        elif input_type == 'ont':
            fastq_input = FastqInput('nanopore', se=[ToolIOFile(input_path)], is_pe=False)
            SnakemakeUtils.dump_object(fastq_input.to_fq_dict(), self._working_dir / 'fq_dict.io')
        else:
            fasta_path = self._working_dir / 'assembly' / 'filtering'
            fasta_path.mkdir(parents=True)
            SnakemakeUtils.dump_object([ToolIOFile(input_path)], fasta_path / 'fasta.io')
        # Config file
        config_data = {
            'contamination_check': {
                'db': str(db),
                'expected_species': expected_species,
                'level_of_depth': level_of_depth},
            'sample_name': sample_name,
            'working_dir': str(self._working_dir),
            'input_type': input_type
        }
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._working_dir)

        # Output files
        output_files = {
            'HTML': self._working_dir / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT,
            'TSV': self._working_dir / contamination_check_kraken.OUTPUT_CONTAMINATION_SUMMARY,
            'INFORMS': self._working_dir / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_INFORMS,
            'INFORMS_KRAKEN': self._working_dir / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS,
        }
        input_format = 'fastq_pe' if input_type == 'illumina' else 'fastq_se' if input_type == 'ont' else 'fasta'
        for k, p in output_files.items():
            if '{input_format}' not in str(p):
                continue
            output_files[k] = Path(str(p).format(input_format=input_format))

        # Run Snakemake
        SnakePipelineUtils.run_snakemake(
            contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN, config_file, list(output_files.values()),
            self._working_dir, threads)
        self.__set_output(output_files)

    def __set_output(self, output_files: Dict[str, Path]) -> None:
        """
        Runs the Snakemake workflow.
        :param output_files: Output files dictionary
        :return: None
        """
        log_file_path = self._working_dir / 'camel.log'
        self._output = Kraken2Output(
            report_section=SnakemakeUtils.load_object(output_files['HTML'])[0].value,
            tsv_summary=output_files['TSV'],
            informs=SnakemakeUtils.load_object(output_files['INFORMS']),
            informs_commands=[SnakemakeUtils.load_object(output_files['INFORMS_KRAKEN'])],
            log_file=log_file_path if log_file_path.exists() else None
        )

    @property
    def output(self) -> Kraken2Output:
        """
        Returns the output of the assembly workflow.
        :return: Assembly output
        """
        return self._output
