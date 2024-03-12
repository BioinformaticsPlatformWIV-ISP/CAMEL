#!/usr/bin/env python
import argparse
import tempfile
from pathlib import Path
from typing import Optional, Dict, List, Sequence

import yaml

from camel.app.camel import Camel
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.loggers import logger
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.pipelines.mycobacterium.bamaddcustomtag import BAMAddCustomTag
from camel.resources.snakefile import variant_calling
from camel.scripts.mycobacteriumpipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainMycobacteriumPipeline(ReportPipeline):
    """
    Main class to run the Mycobacterium pipeline.
    """

    CUSTOM_ANALYSES = [
        'kraken2', 'ncbi_16s', 'csb_rd', 'hsp65', '51snp', 'snpit', 'spoligotyping', 'snp_lineage', 'amr', 'mlst',
        'cgmlst', 'confindr']

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Mycobacterium pipeline', '1.2', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return '<i>Mycobacterium</i> pipeline'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        self._validate_input_files()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)
        self._export_assembly()
        self._export_bam()

    def _export_bam(self) -> None:
        """
        Exports the BAM output to the specified output location (optional).
        :return: None
        """
        if self._args.output_bam is None:
            logger.debug(f'Not exporting BAM output')
            return
        path_io = self._args.working_dir / variant_calling.OUTPUT_VARIANT_CALLING_BAM
        bam_input = SnakemakeUtils.load_object(path_io)

        # Add custom tag with the sample name for PACU
        with tempfile.TemporaryDirectory(prefix='camel_', dir=Camel.get_instance().config['temp_dir']) as dir_:
            add_tag = BAMAddCustomTag(Camel.get_instance())
            add_tag.add_input_files({'BAM': bam_input})
            add_tag.update_parameters(output=str(self._args.output_bam), name='PACU_name', value=self.sample_name)
            add_tag.run(Path(str(dir_)))
        logger.info(f'Output BAM file copied to: {self._args.output_bam}')

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments (optional)
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        parser.add_argument(
            '--output-bam', type=Path, help='Output path for the mapping to the reference genome (BAM)')
        for analysis_key in MainMycobacteriumPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        return parser.parse_args(args)

    def __construct_config_file(self, input_files: Dict[str, List[Dict[str, str]]]) -> str:
        """
        Constructs the configuration file.
        :param input_files: Dictionary with the input files (keys can be FASTQ_PE, FASTQ_SE).
        :return: Configuration file
        """
        config_data = self.get_template_data(input_files)
        config_data['analyses'] = [key for key in MainMycobacteriumPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with open(CONFIG_DATA) as handle_in:
            config_data.update(yaml.load(handle_in.read().format(
                qc_typing_scheme='cgmlst' if self._args.cgmlst else 'mlst',
                export_fastq='true' if self._args.report_include_fastq else 'false',
                export_bam='true' if self._args.report_include_bam else 'false',
                coverage_max=self._args.cov_max
            ), Loader=yaml.SafeLoader))
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainMycobacteriumPipeline()
    main.run()
