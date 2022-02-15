import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Any

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
from camel.resources.snakefile import variant_filtering, variant_calling


@dataclass
class VariantFilteringOutput:
    """
    This dataclass holds the output of the variant filtering workflow.
    """
    vcf_filtered: ToolIOFile
    stats: Dict
    nb_of_variants: int
    informs: List[Dict[str, Any]]


class VariantFilteringWrapper(object):
    """
    This class is used as a wrapper class around the variant filtering Snakemake workflow.
    """

    def __init__(self, working_dir: Path) -> None:
        """
        Initializes the variant calling wrapper.
        :param working_dir: Working directory
        """
        self._working_dir = working_dir
        self._output = None

    @property
    def output(self) -> VariantFilteringOutput:
        """
        Returns the filtered VCF file
        :return: VCF file path
        """
        return self._output

    def __convert_to_vcf_gz(self, vcf_file: Path) -> Path:
        """
        Converts an input VCF file to an indexed VCF_GZ file.
        :param vcf_file: Input VCF file
        :return: Indexed VCF_GZ file
        """
        c = Camel()
        bcftools_view = BcftoolsView(c)
        bcftools_view.add_input_files({'VCF': [ToolIOFile(vcf_file)]})
        input_dir = Path(self._working_dir, 'input')
        if not input_dir.is_dir():
            input_dir.mkdir(parents=True)
        bcftools_view.update_parameters(compress_output=True, output_format='VCF')
        bcftools_view.run(input_dir)
        return bcftools_view.tool_outputs['VCF_GZ'][0].path

    def __create_input(self, vcf_gz_file: Path, bam_file: Path) -> None:
        """
        Creates the input files for the workflow.
        :param vcf_gz_file: Input VCF GZ file
        :param bam_file: Input BAM file
        :return: None
        """
        for path, destination in [(vcf_gz_file, variant_calling.OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ),
                                  (bam_file, variant_calling.OUTPUT_VARIANT_CALLING_BAM)]:
            target_file = Path(self._working_dir, destination)
            if not target_file.parent.exists():
                target_file.parent.mkdir(parents=True)
            SnakemakeUtils.dump_object([ToolIOFile(path)] if path is not None else [], target_file)

    def run_workflow(
            self, sample_name: str, vcf_file: Path, bam_file: Path, filtering_options: Dict, cores: int = 8) -> None:
        """
        Runs the variant calling workflow.
        :param sample_name: Sample name
        :param vcf_file: Input VCF file
        :param bam_file: Input BAM file
        :param cores: Number of cores
        :param filtering_options: Dict
        :return: None
        """
        if not self._working_dir.exists():
            self._working_dir.mkdir(parents=True)
        vcf_gz_file = self.__convert_to_vcf_gz(vcf_file)
        self.__create_input(vcf_gz_file, bam_file)

        # Create config
        config_path = SnakePipelineUtils.generate_config_file({
            'working_dir': str(self._working_dir),
            'variant_filtering': filtering_options,
            'sample_name': sample_name
        }, self._working_dir)

        # Execute Snakemake
        output_files = {
            'VCF': self._working_dir / variant_filtering.OUTPUT_VARIANT_FILTERING_VCF,
            'STATS': self._working_dir / variant_filtering.OUTPUT_VARIANT_FILTERING_STATS,
            'INFORMS': self._working_dir / variant_filtering.OUTPUT_VARIANT_FILTERING_INFORMS_ALL
        }
        SnakePipelineUtils.run_snakemake(
            variant_filtering.SNAKEFILE_VARIANT_FILTERING, config_path, list(output_files.values()), self._working_dir,
            cores)
        self.__set_output(output_files)

    def __set_output(self, output_files: Dict[str, Path]) -> None:
        """
        Collects the output of the workflow.
        :param output_files: Output files
        :return: None
        """
        json_file = SnakemakeUtils.load_object(output_files['STATS'])[0].path
        with open(json_file) as handle:
            stats = json.load(handle)
        self._output = VariantFilteringOutput(
            vcf_filtered=SnakemakeUtils.load_object(output_files['VCF'])[0],
            stats=stats,
            nb_of_variants=stats['zscore']['variants_out'],
            informs=SnakemakeUtils.load_object(output_files['INFORMS'])
        )
