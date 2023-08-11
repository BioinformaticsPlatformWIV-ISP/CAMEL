import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import quality_checks, contamination_check_kraken, medaka_polishing, assembly_spades, \
    variant_calling, trimming_illumina, trimming_ont


@dataclass
class QCOutput:
    report_section: HtmlReport
    tsv_summary: Path
    # informs: Path


class QCWrapper(object):
    """
    This class is used as a wrapper class around the quality checks workflow.
    """

    def __init__(self, working_dir: Path) -> None:
        """
        Initializes the quality checks wrapper.
        :param working_dir: Working directory
        :return: None
        """
        self._working_dir = working_dir
        self._output = None

    def run_workflow(
            self, trimming_input: List[Path], read_type: str, config_parameters: Dict, kraken_input=None,
            mapping_rate_input=None, depth_informs=None, typing_informs=None, nanopore_specific=None,
            threads: int = 8) -> None:
        """
        Runs the QC workflow.
        :param nanopore_specific: Reads and nanoplot informs output
        :param trimming_input: Fastqc or nanoplot output
        :param typing_informs: IO file of the typing informs
        :param depth_informs: Depth informs of spades or medaka
        :param mapping_rate_input: Mapping rate of medaka or spades
        :param kraken_input: Kraken2 informs
        :param read_type: Sequencing technology
        :param config_parameters: Config file parameters
        :param threads: Number of threads
        :return: None
        """
        if kraken_input is None:
            kraken_input = []
        if typing_informs is None:
            typing_informs = []
        if depth_informs is None:
            depth_informs = []
        if mapping_rate_input is None:
            mapping_rate_input = []
        if nanopore_specific is None:
            nanopore_specific = []
        output_directory = self._working_dir / 'qc_{}'.format(read_type)
        config_data = {
            'working_dir': str(output_directory),
            'read_type': read_type,
            'analyses': config_parameters['analyses'],
            'quality_checks': {'coverage_mode': config_parameters.get('coverage_mode', 'assembly'),
                               'typing_scheme': config_parameters['quality_checks'].get('typing_scheme', None),
                               'expected_gc_content': config_parameters['quality_checks']['expected_gc_content'],
                               'disabled_checks': config_parameters['quality_checks'].get('disabled_checks', [])}
        }

        # for debug - needs to think about it
        config_data['quality_checks']['disabled_checks'].append('typing')

        config_file = SnakePipelineUtils.generate_config_file(config_data, output_directory)

        kraken_pickle = Path(output_directory) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_INFORMS

        mapping_rate_informs_nanopore = Path(output_directory) / medaka_polishing.OUTPUT_ASSEMBLY_MAPPING_RATE_INFORMS
        mapping_rate_informs_illumina = Path(output_directory) / assembly_spades.OUTPUT_ASSEMBLY_MAPPING_INFORMS
        mapping_rate_informs_illumina_variant = Path(
            output_directory) / variant_calling.OUTPUT_VARIANT_CALLING_MAPPING_INFORMS

        depth_informs_nanopore = Path(output_directory) / medaka_polishing.OUTPUT_ASSEMBLY_DEPTH_INFORMS
        depth_informs_illumina = Path(output_directory) / assembly_spades.OUTPUT_ASSEMBLY_DEPTH_INFORMS
        depth_informs_illumina_variant = Path(output_directory) / variant_calling.OUTPUT_VARIANT_CALLING_DEPTH_INFORMS

        typing_informs_path = Path(output_directory) / 'typing' / '{typing_wildcard}' / 'stats' / 'informs.io'

        trimming_illumina_pre = Path(output_directory) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_FASTQC_TXT_PRE
        trimming_illumina_post = Path(output_directory) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_FASTQC_TXT_POST

        trimming_nanopore_pre = Path(output_directory) / trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_TXT_PRE
        trimming_nanopore_post = Path(output_directory) / trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_TXT_POST

        nanopore_reads = Path(output_directory) / trimming_ont.OUTPUT_TRIMMING_ONT_READS
        nanopore_nanoplot_informs = Path(output_directory) / trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_INFORMS_POST

        if read_type == 'illumina':
            mode = config_data['quality_checks']['coverage_mode']
            if mode == 'assembly':
                mapping_rate_informs_illumina.parent.mkdir(exist_ok=True, parents=True)
                shutil.copyfile(mapping_rate_input, mapping_rate_informs_illumina)
                depth_informs_illumina.parent.mkdir(exist_ok=True, parents=True)
                shutil.copyfile(depth_informs, depth_informs_illumina)
            else:
                mapping_rate_informs_illumina_variant.parent.mkdir(exist_ok=True, parents=True)
                shutil.copyfile(mapping_rate_input, mapping_rate_informs_illumina_variant)
                depth_informs_illumina_variant.parent.mkdir(exist_ok=True, parents=True)
                shutil.copyfile(depth_informs, depth_informs_illumina)
            trimming_illumina_pre.parent.mkdir(exist_ok=True, parents=True)
            trimming_illumina_post.parent.mkdir(exist_ok=True, parents=True)
            shutil.copyfile(trimming_input[0], trimming_illumina_pre)
            shutil.copyfile(trimming_input[1], trimming_illumina_post)
        elif read_type == 'nanopore':
            mapping_rate_informs_nanopore.parent.mkdir(exist_ok=True, parents=True)
            depth_informs_nanopore.parent.mkdir(exist_ok=True, parents=True)
            trimming_nanopore_pre.parent.mkdir(exist_ok=True, parents=True)
            trimming_nanopore_post.parent.mkdir(exist_ok=True, parents=True)
            nanopore_reads.parent.mkdir(exist_ok=True, parents=True)
            nanopore_nanoplot_informs.parent.mkdir(exist_ok=True, parents=True)
            shutil.copyfile(mapping_rate_input, mapping_rate_informs_nanopore)
            shutil.copyfile(depth_informs, depth_informs_nanopore)
            shutil.copyfile(trimming_input[0], trimming_nanopore_pre)
            shutil.copyfile(trimming_input[1], trimming_nanopore_post)
            shutil.copyfile(nanopore_specific[0], nanopore_reads)
            shutil.copyfile(nanopore_specific[1], nanopore_nanoplot_informs)
        if 'kraken' in config_data['analyses']:
            kraken_pickle.parent.mkdir(exist_ok=True, parents=True)
            shutil.copyfile(kraken_input, kraken_pickle)

        # for mlst_key in [key for key in config_data['analyses'] if 'mlst' in key]:
        #     typing_under_study = typing_informs_path.format(typing_wildcard=mlst_key)
        #     typing_under_study.parent.mkdir(exist_ok=True, parents=True)
        #     shutil.copyfile(typing_informs[mlst_key], typing_under_study)

        output_files = {
            'HTML': output_directory / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT,
            'TSV': output_directory / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
            # 'INFORMS': self._working_dir / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT_JSON,
        }
        SnakePipelineUtils.run_snakemake(
            quality_checks.SNAKEFILE_QUALITY_CHECKS, config_file, list(output_files.values()),
            output_directory, threads)
        self.__set_output(output_files)

    def __set_output(self, output_files: Dict[str, Path]) -> None:
        """
        Runs the Snakemake workflow.
        :param output_files: Output files dictionary
        :return: None
        """
        self._output = QCOutput(
            report_section=SnakemakeUtils.load_object(output_files['HTML'])[0].value,
            tsv_summary=output_files['TSV'],
            # informs=output_files['INFORMS']
        )

    @property
    def output(self) -> QCOutput:
        """
        Returns the output of the assembly workflow.
        :return: Assembly output
        """
        return self._output
