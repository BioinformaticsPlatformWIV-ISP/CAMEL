#!/usr/bin/env python
import argparse
import json
import logging
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.components.workflows.genedetectionwrapper import GeneDetectionWrapper, GeneDetectionOutput


class MainGeneDetection(object):
    """
    This class is used to run the gene detection tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments (optional)
        """
        self._args = MainGeneDetection.parse_arguments(args)
        self._sample_name = MainScriptHelper.determine_sample_name(self._args)
        self._helper = MainScriptHelper(self._args.working_dir, self._sample_name)
        self._report = None

    @staticmethod
    def parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        MainScriptHelper.add_common_arguments(argument_parser)
        MainScriptHelper.add_assembly_arguments(argument_parser)
        MainScriptHelper.add_input_files_arguments(argument_parser)
        group_db = argument_parser.add_mutually_exclusive_group(required=True)
        group_db.add_argument('--database-dir', type=str)
        group_db.add_argument('--database-html', type=str)
        argument_parser.add_argument('--detection-method', type=str, choices=['blast', 'srst2', 'kma'], default='blast')

        # BLAST specific parameters
        argument_parser.add_argument('--blast-min-percent-identity', type=int, default=90)
        argument_parser.add_argument('--blast-min-percent-coverage', type=int, default=60)
        argument_parser.add_argument('--blast-task', type=str, choices=['blastn', 'megablast'], default='megablast')

        # SRST2 specific parameters
        argument_parser.add_argument('--srst2-min-cov', type=int, default=90)
        argument_parser.add_argument('--srst2-max-div', type=int, default=10)
        argument_parser.add_argument('--srst2-max-unaligned-overlap', type=int, default=100)
        argument_parser.add_argument('--srst2-max-mismatch', type=int, default=10)

        # KMA specific parameters
        argument_parser.add_argument('--kma-min-percent-identity', type=int, default=90)
        argument_parser.add_argument('--kma-min-percent-coverage', type=int, default=60)
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        self._report = self._helper.init_report(
            self._args.output_html, self._args.output_dir, 'Gene detection report',
            f'Gene detection {self._args.detection_method}')
        self.__add_analysis_info_section()
        input_files = self._helper.symlink_input_files(self._args.fasta, self._args.fastq_pe)
        db_data = self.__get_db_metadata()

        wrapper = GeneDetectionWrapper(self._args.working_dir)
        if self._args.detection_method == 'blast':
            fasta_file = self._helper.get_blast_input(input_files, self._args, self._report).path
            wrapper.run_workflow_blast(fasta_file, self._sample_name, db_data, self._args.threads)
        elif self._args.detection_method == 'srst2':
            fq_files = [f.path for f in self._helper.get_srst2_input(input_files, self._args, self._report)]
            wrapper.run_workflow_srst2(fq_files, self._sample_name, db_data, self._args.threads)
        else:
            fq_files = [f.path for f in self._helper.get_srst2_input(input_files, self._args, self._report)]
            wrapper.run_workflow_kma(fq_files, self._sample_name, db_data, self._args.threads)
        self.__export_output(wrapper.output)

    def __add_analysis_info_section(self) -> None:
        """
        Adds the report section with the analysis info
        :return: None
        """
        input_files_str = self._helper.determine_input_files(self._args)
        self._helper.export_analysis_info_section(self._report, input_files_str)

    def __get_db_metadata(self) -> Dict[str, Any]:
        """
        Returns the database information dictionary.
        :return: Database information dictionary
        """
        # Get database path
        if self._args.database_dir is not None:
            db_path = Path(self._args.database_dir)
        else:
            db_path = Path(f"{'.'.join(self._args.database_html.split('.')[:-1])}_files")
        config_data = {'path': db_path}

        # Add specific options
        if self._args.detection_method == 'blast':
            config_data.update({'params': {'blastn': {
                'min_percent_identity': self._args.blast_min_percent_identity,
                'min_coverage': self._args.blast_min_percent_coverage,
                'task': self._args.blast_task
            }}})
        elif self._args.detection_method == 'srst2':
            config_data.update({'params': {'srst2': {
                'min_coverage': self._args.srst2_min_cov,
                'max_divergence': self._args.srst2_max_div,
                'max_unaligned_overlap': self._args.srst2_max_unaligned_overlap,
                'max_mismatch': self._args.srst2_max_mismatch
            }}})
        elif self._args.detection_method == 'kma':
            config_data.update({'params': {'kma': {
                'min_percent_identity': self._args.kma_min_percent_identity,
                'min_coverage': self._args.kma_min_percent_coverage,
            }}})

        # Add extra column
        with (db_path / 'db_metadata.txt').open() as handle:
            db_metadata = json.load(handle)
            if 'extra_column' in db_metadata:
                config_data['metadata'] = db_metadata['extra_column']
        return config_data

    def __export_output(self, output: GeneDetectionOutput) -> None:
        """
        Exports the output of the workflow.
        :param output: Output
        :return: None
        """
        self._helper.logs['gene_detection'] = str(output.log_file)
        self._helper.informs.append(output.informs)
        self._helper.export_output_and_commands_section(self._report, output.report_section)


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main = MainGeneDetection()
    main.run()
