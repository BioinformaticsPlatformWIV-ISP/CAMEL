import json
from dataclasses import dataclass
from typing import Dict

import os
import yaml

from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
from camel.resources.snakefile import SNAKEFILE_VARIANT_FILTERING
from camel.resources.snakefile.variant_calling import OUTPUT_VARIANT_CALLING_BAM, \
    OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ
from camel.resources.snakefile.variant_filtering import OUTPUT_VARIANT_FILTERING_VCF, OUTPUT_VARIANT_FILTERING_STATS


class VariantFilteringWrapper(object):
    """
    This class is used as a wrapper class around the variant filtering Snakemake workflow.
    """

    @dataclass
    class VariantFilteringOutput:
        vcf_filtered: ToolIOFile
        stats: Dict
        nb_of_variants: int

    def __init__(self, working_dir: str) -> None:
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

    def __convert_to_vcf_gz(self, vcf_file) -> str:
        """
        Converts an input VCF file to an indexed VCF_GZ file.
        :param vcf_file: Input VCF file
        :return: Indexed VCF_GZ file
        """
        c = Camel()
        bcftools_view = BcftoolsView(c)
        bcftools_view.add_input_files({'VCF': [ToolIOFile(vcf_file)]})
        input_dir = os.path.join(self._working_dir, 'input')
        if not os.path.isdir(input_dir):
            os.makedirs(input_dir)
        bcftools_view.update_parameters(compress_output=True, output_format='VCF')
        bcftools_view.run(input_dir)
        return bcftools_view.tool_outputs['VCF_GZ'][0].path

    def __create_input(self, vcf_gz_file: str, bam_file: str) -> None:
        """
        Creates the input files for the workflow.
        :param vcf_gz_file: Input VCF GZ file
        :param bam_file: Input BAM file
        :return: None
        """
        for path, destination in [(vcf_gz_file, OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ),
                                  (bam_file, OUTPUT_VARIANT_CALLING_BAM)]:
            target_dir = os.path.dirname(os.path.join(self._working_dir, destination))
            if not os.path.isdir(target_dir):
                os.makedirs(target_dir)
            SnakemakeUtils.dump_object([ToolIOFile(path)] if path is not None else [],
                                       os.path.join(self._working_dir, destination))

    def run_workflow(self, vcf_file: str, bam_file: str, filtering_options: Dict, cores: int = 8) -> None:
        """
        Runs the variant calling workflow.
        :param vcf_file: Input VCF file
        :param bam_file: Input BAM file
        :param cores: Number of cores
        :param filtering_options: Dict
        :return: None
        """
        if not os.path.isdir(self._working_dir):
            os.makedirs(self._working_dir)
        vcf_gz_file = self.__convert_to_vcf_gz(vcf_file)
        self.__create_input(vcf_gz_file, bam_file)

        # Create config
        config_path = os.path.join(self._working_dir, 'config_variant_filtering.yaml')
        with open(config_path, 'w') as handle:
            yaml.dump({
                'working_dir': self._working_dir,
                'variant_filtering': filtering_options
            }, handle)

        # Execute Snakemake
        output_files = {
            'VCF': os.path.join(self._working_dir, OUTPUT_VARIANT_FILTERING_VCF),
            'STATS': os.path.join(self._working_dir, OUTPUT_VARIANT_FILTERING_STATS)
        }
        command = Command('snakemake --snakefile {} --configfile {} {} --cores {}'.format(
            SNAKEFILE_VARIANT_FILTERING, config_path, ' '.join(output_files.values()), cores
        ))
        command.run_command(self._working_dir)
        self.__set_output(output_files)

    def __set_output(self, output_files: Dict[str, str]) -> None:
        """
        Collects the output of the workflow.
        :param output_files: Output files
        :return: None
        """
        json_file = SnakemakeUtils.load_object(output_files['STATS'])[0].path
        with open(json_file) as handle:
            stats = json.load(handle)
        self._output = VariantFilteringWrapper.VariantFilteringOutput(
            vcf_filtered=SnakemakeUtils.load_object(output_files['VCF'])[0],
            stats=stats,
            nb_of_variants=stats['zscore']['variants_out']
        )
