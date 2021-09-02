#!/usr/bin/env python
import argparse
import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Sequence

import yaml

from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.stecpipeline import CONFIG_DATA
from camel.scripts.stecpipeline import SNAKEFILE_MAIN


class MainSTECPipeline(ReportPipeline):
    """
    Main class to run the STEC pipeline.
    """

    CUSTOM_ANALYSES = ['kraken', 'resfinder', 'argannot', 'card', 'ncbi_amr', 'mlst_pasteur', 'mlst_warwick', 'cgmlst',
                       'pointfinder', 'plasmidfinder', 'serotype', 'virulencefinder', 'innuendo_cgmlst']

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('STEC pipeline', '1.0', SNAKEFILE_MAIN, args)

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)

    def __construct_config_file(self, input_files: List[Dict[str, str]]) -> str:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        key = 'fastq_pe' if (self._args.fastq_pe is not None) else 'fastq_se'
        config_data = self.get_template_data(key, input_files)
        with open(CONFIG_DATA) as handle_in:
            config_data.update(yaml.load(handle_in.read(), Loader=yaml.SafeLoader))
        config_data['analyses'] = [key for key in MainSTECPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        config_data['read_type'] = self._args.read_type
        config_data['quality_checks']['typing_scheme'] = 'cgmlst' if self._args.cgmlst else 'mlst_warwick'
        config_data['read_trimming']['export_fastq'] = 'true' if self._args.report_include_fastq else 'false'
        if self._args.library is not None:
            config_data['read_trimming']['adapter'] = self._args.library
        config_data['variant_calling']['report_include_bam'] = 'true' if self._args.report_include_bam else 'false'
        detection_method_cgmlst = {
            'blast': 'blast', 'srst2': 'blast', 'kma': 'kma'}.get(self._args.detection_method)
        config_data['sequence_typing']['cgmlst']['detection_method'] = detection_method_cgmlst
        config_data['sequence_typing']['innuendo_cgmlst']['detection_method'] = detection_method_cgmlst
        if self._args.read_type == 'iontorrent':
            config_data['assembly']['spades']['iontorrent'] = None
        return SnakePipelineUtils.generate_config_file(config_data, Path(self._args.working_dir))

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        parser.add_argument('--fastq-se', help="Input SE FASTQ file")
        parser.add_argument('--fastq-se-name', help="Input SE FASTQ file name")
        parser.add_argument(
            '--read-type', help="Type of reads.", choices=['illumina', 'iontorrent'], default='illumina')
        for analysis_key in MainSTECPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        return parser.parse_args(args)

    def _get_fastq_input_links(self) -> List[List[Tuple[str, str]]]:
        """
        Returns the links to the input FASTQ files.
        :return: Links
        """
        links = []
        if self._args.fastq_se is not None:
            gzipped = FileSystemHelper.is_gzipped(self._args.fastq_se)
            links.append([self._args.fastq_se, f"{self.sample_name}.fastq{'.gz' if gzipped else ''}"])
        else:
            for read_nb, path in enumerate(self._args.fastq_pe, start=1):
                gzipped = FileSystemHelper.is_gzipped(path)
                links.append([path, f"{self.sample_name}_{read_nb}.fastq{'.gz' if gzipped else ''}"])
        return links

    @property
    def sample_name(self) -> str:
        """
        Returns the sample name.
        :return: Sample name
        """
        if self._args.fastq_pe is not None:
            return super().sample_name
        else:
            name = self._args.fastq_se_name if (self._args.fastq_se_name is not None) else self._args.fastq_se
            return FastqUtils.get_sample_name(name, FastqUtils.PATTERN_FQ_SE)


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main = MainSTECPipeline()
    main.run()
