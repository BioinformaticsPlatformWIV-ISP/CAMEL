#!/usr/bin/env python
import argparse
import shutil
from typing import Optional, List, Dict, Sequence

import yaml

from camel.app.camel import Camel
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.loggers import logger
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import assembly_spades
from camel.scripts.listeriapipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainListeriaPipeline(ReportPipeline):
    """
    Main class to run the Listeria pipeline.
    """

    CUSTOM_ANALYSES = [
        'kraken', 'confindr', 'mlst', 'cgmlst', 'species_confirmation', 'resfinder', 'argannot', 'card', 'ncbi_amr',
        'virulencefinder', 'vfdb_core', 'plasmidfinder', 'typing_amr', 'typing_virulence', 'pcr_serogroup']

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Listeria pipeline', '1.3', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return '<i>Listeria</i> pipeline'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        self._validate_fastq_input()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)

        # TODO: move to base class
        if self._args.output_fasta is not None:
            path_io = self._args.working_dir / assembly_spades.OUTPUT_ASSEMBLY_FASTA
            path_fasta = SnakemakeUtils.load_object(path_io)[0].path
            shutil.copyfile(path_fasta, self._args.output_fasta)
            logger.info(f'Output FASTA file copied to: {self._args.output_fasta}')

    def __construct_config_file(self, input_files: List[Dict[str, str]]) -> str:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        config_data = self.get_template_data('fastq_pe', input_files)
        config_data['analyses'] = [key for key in MainListeriaPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with CONFIG_DATA.open() as handle_in:
            config_data.update(yaml.load(handle_in.read().format(
                qc_typing_scheme='cgmlst' if self._args.cgmlst else 'mlst',
                export_fastq='true' if self._args.report_include_fastq else 'false',
                export_bam='true' if self._args.report_include_bam else 'false',
                coverage_max=self._args.cov_max
            ), Loader=yaml.SafeLoader))
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        for analysis_key in MainListeriaPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        return parser.parse_args(args)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainListeriaPipeline()
    main.run()
