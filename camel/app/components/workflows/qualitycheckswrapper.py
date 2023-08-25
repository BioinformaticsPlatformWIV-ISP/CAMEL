import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import quality_checks, contamination_check_kraken, medaka_polishing, assembly_spades, \
    variant_calling, trimming_illumina, trimming_ont


@dataclass
class QCOutput:
    report_section: HtmlReportSection
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
        if not self._working_dir.exists():
            self._working_dir.mkdir(parents=True)
        config_data = {
            'working_dir': str(self._working_dir),
            'read_type': read_type,
            'analyses': config_parameters['analyses'],
            'quality_checks': {'coverage_mode': config_parameters.get('coverage_mode', 'assembly'),
                               'typing_scheme': config_parameters['quality_checks'].get('typing_scheme', None),
                               'expected_gc_content': config_parameters['quality_checks']['expected_gc_content'],
                               'disabled_checks': config_parameters['quality_checks'].get('disabled_checks', [])}
        }

        # for debug - needs to think about it
        # config_data['quality_checks']['disabled_checks'].append('typing')

        kraken_pickle = Path(self._working_dir) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_INFORMS

        mapping_rate_informs_nanopore = Path(self._working_dir) / medaka_polishing.OUTPUT_ASSEMBLY_MAPPING_RATE_INFORMS
        mapping_rate_informs_illumina = Path(self._working_dir) / assembly_spades.OUTPUT_ASSEMBLY_MAPPING_INFORMS
        mapping_rate_informs_illumina_variant = Path(
            self._working_dir) / variant_calling.OUTPUT_VARIANT_CALLING_MAPPING_INFORMS

        depth_informs_nanopore = Path(self._working_dir) / medaka_polishing.OUTPUT_ASSEMBLY_DEPTH_INFORMS
        depth_informs_illumina = Path(self._working_dir) / assembly_spades.OUTPUT_ASSEMBLY_DEPTH_INFORMS
        depth_informs_illumina_variant = Path(self._working_dir) / variant_calling.OUTPUT_VARIANT_CALLING_DEPTH_INFORMS

        trimming_illumina_pre = Path(self._working_dir) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_FASTQC_TXT_PRE
        trimming_illumina_post = Path(self._working_dir) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_FASTQC_TXT_POST

        trimming_nanopore_pre = Path(self._working_dir) / trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_TXT_PRE
        trimming_nanopore_post = Path(self._working_dir) / trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_TXT_POST

        nanopore_reads = Path(self._working_dir) / trimming_ont.OUTPUT_TRIMMING_ONT_READS
        nanopore_nanoplot_informs = Path(self._working_dir) / trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_INFORMS_POST

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

        if 'cgmlst' in config_data['quality_checks']['typing_scheme'] and typing_informs is not None:
            typing_informs_path = str(Path(self._working_dir) / 'typing' / '{typing_wildcard}' / 'stats' /
                                      'informs.io').format(typing_wildcard=config_data['quality_checks']['typing_scheme'])
            Path(typing_informs_path).parent.mkdir(exist_ok=True, parents=True)
            shutil.copyfile(typing_informs, typing_informs_path)
        else:
            config_data['quality_checks']['disabled_checks'] = 'typing'

        config_file = SnakePipelineUtils.generate_config_file(config_data, self._working_dir)

        output_files = {
            'HTML': self._working_dir / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT,
            'TSV': self._working_dir / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
            # 'INFORMS': self._working_dir / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT_JSON,
        }
        SnakePipelineUtils.run_snakemake(
            quality_checks.SNAKEFILE_QUALITY_CHECKS, config_file, list(output_files.values()),
            self._working_dir, threads)
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
