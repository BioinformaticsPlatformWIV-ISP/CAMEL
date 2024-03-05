#!/usr/bin/env python
import argparse
from typing import Tuple, Any, List, Dict, Optional, Sequence

import yaml

from camel.app.camel import Camel
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.shigellapipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainShigellaPipeline(ReportPipeline):
    """
    Main class to run the Shigella pipeline.
    """

    CUSTOM_ANALYSES = [
        'kraken', 'confindr', 'resfinder', 'argannot', 'card', 'ncbi_amr', 'mlst_pasteur', 'mlst_warwick', 'cgmlst',
        'pointfinder', 'plasmidfinder', 'virulencefinder', 'identification']

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Shigella pipeline', '0.6', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return '<i>Shigella</i> pipeline'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        self._validate_input_files()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)

    def __construct_config_file(self, input_files: List[Dict[str, str]]) -> str:
        """
        Constructs the configuration file to run Snakemake.
        :return: Config file path
        """
        config_data = self.get_template_data('fastq_pe', input_files)
        config_data['analyses'] = [key for key in MainShigellaPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with open(CONFIG_DATA) as handle_in:
            config_data.update(yaml.safe_load(handle_in.read().format(
                export_fastq='true' if self._args.report_include_fastq else "'false'",
                export_bam='true' if self._args.report_include_bam else 'false',
                variant_filtering={},
                coverage_max=self._args.cov_max
            )))
        config_data['quality_checks']['typing_scheme'] = 'cgmlst' if self._args.cgmlst else 'mlst_warwick'
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    def __create_fastq_input_dict(self) -> Tuple[str, List[Dict[str, Any]]]:
        """
        Creates the input dictionary with the FASTQ files.
        :return: Input key, input dictionary
        """
        fq_files = SnakePipelineUtils.symlink_input_files(
            self._args.working_dir / 'input', self._args.fastq_pe, self._args.fastq_pe_names, True)
        return 'fastq_pe', [{'name': name, 'path': path} for name, path in zip(
            self._args.fastq_pe_names, fq_files)]

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        for analysis_key in MainShigellaPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        return parser.parse_args(args)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainShigellaPipeline()
    main.run()
