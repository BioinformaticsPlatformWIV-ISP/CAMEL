#!/usr/bin/env python
import argparse
import json
import logging
from typing import Any, Dict, List, Optional

import os

from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.components.workflows.genedetectionwrapper import GeneDetectionWrapper, GeneDetectionOutput
from camel.app.io.tooliofile import ToolIOFile


class MainGeneDetection(object):
    """
    This class is used to run the gene detection tool.
    """

    def __init__(self, args: Optional[argparse.Namespace] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments (optional)
        """
        self._args = args if args is not None else MainGeneDetection.parse_arguments()
        self._sample_name = MainScriptHelper.determine_sample_name(self._args)
        self._helper = MainScriptHelper(self._args.working_dir, self._sample_name)
        self._report = None

    @staticmethod
    def parse_arguments() -> argparse.Namespace:
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
        argument_parser.add_argument('--detection-method', type=str, choices=['blast', 'srst2'], default='blast')
        argument_parser.add_argument('--report-include-fastq', action='store_true')

        # BLAST specific parameters
        argument_parser.add_argument('--blast-min-percent-identity', type=int, default=90)
        argument_parser.add_argument('--blast-min-percent-coverage', type=int, default=60)

        # SRST2 specific parameters
        argument_parser.add_argument('--srst2-min-cov', type=int, default=90)
        argument_parser.add_argument('--srst2-max-div', type=int, default=10)
        argument_parser.add_argument('--srst2-max-unaligned-overlap', type=int, default=100)
        argument_parser.add_argument('--srst2-max-mismatch', type=int, default=10)
        return argument_parser.parse_args()

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
        if self._args.detection_method == 'blast':
            fasta_file = self._helper.get_blast_input(input_files, self._args, self._report)
            output = self.__run_gene_detection_blast(fasta_file, db_data)
        else:
            fastq_files = self._helper.get_srst2_input(input_files, self._args, self._report)
            output = self.__run_gene_detection_srst2(fastq_files, db_data)
        self.__export_output(output)

    def __add_analysis_info_section(self) -> None:
        """
        Adds the report section with the analysis info
        :return: None
        """
        if self._args.detection_method == 'blast':
            data = [
                ['% identity threshold:', f'{self._args.blast_min_percent_identity}%'],
                ['% query covered threshold:', f'{self._args.blast_min_percent_coverage}%']
            ]
        elif self._args.detection_method == 'srst2':
            data = [
                ['Min. % coverage threshold:', f'{self._args.srst2_min_cov}%'],
                ['Max. % divergence threshold:', f'{self._args.srst2_max_div}%']
            ]
        else:
            raise ValueError(f"Invalid detection method: {self._args.detection_method}")
        input_files_str = self._helper.determine_input_files(self._args)
        self._helper.export_analysis_info_section(self._report, input_files_str, data)

    def __get_db_metadata(self) -> Dict[str, Any]:
        """
        Returns the database information dictionary.
        :return: Database information dictionary
        """
        # Get database path
        if self._args.database_dir is not None:
            db_path = self._args.database_dir
        else:
            db_path = f"{'.'.join(self._args.database_html.split('.')[:-1])}_files"
        metadata = {'path': db_path}

        # Add specific options
        if self._args.detection_method == 'blast':
            metadata.update({'blast_filtering_options': {
                'min_percent_identity': self._args.blast_min_percent_identity,
                'min_coverage': self._args.blast_min_percent_coverage
            }})
        elif self._args.detection_method == 'srst2':
            metadata.update({'srst2_options': {
                'min_coverage': self._args.srst2_min_cov,
                'max_divergence': self._args.srst2_max_div,
                'max_unaligned_overlap': self._args.srst2_max_unaligned_overlap,
                'max_mismatch': self._args.srst2_max_mismatch
            }})

        # Add extra column
        with open(os.path.join(db_path, 'db_metadata.txt')) as handle:
            db_metadata = json.load(handle)
            if 'extra_column' in db_metadata:
                metadata['extra_column'] = db_metadata['extra_column']
        return metadata

    def __run_gene_detection_blast(self, fasta_file: ToolIOFile, db_data: Dict[str, Any]) -> GeneDetectionOutput:
        """
        Runs the gene detection workflow.
        :param fasta_file: FASTA file
        :param db_data: Database information dictionary
        :return: None
        """
        wrapper = GeneDetectionWrapper(self._args.working_dir)
        wrapper.run_workflow_blast(fasta_file.path, self._sample_name, db_data, self._args.threads)
        return wrapper.output

    def __run_gene_detection_srst2(self, fastq_pe: List[ToolIOFile], db_data: Dict[str, Any]) -> GeneDetectionOutput:
        """
        Runs the gene detection workflow in srst2 mode.
        :param fastq_pe: Paired end FASTQ input
        :param db_data: Database information dictionary
        :return: None
        """
        wrapper = GeneDetectionWrapper(os.path.join(self._args.working_dir, os.path.basename(db_data['path'])))
        wrapper.run_workflow_srst2([f.path for f in fastq_pe], self._sample_name, db_data, self._args.threads)
        return wrapper.output

    def __export_output(self, output: GeneDetectionOutput) -> None:
        """
        Exports the output of the workflow.
        :param output: Output
        :return: None
        """
        self._helper.logs['gene_detection'] = output.log_file
        self._helper.informs.append(output.informs)
        self._helper.export_output_and_commands_section(self._report, output.report_section)


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main = MainGeneDetection()
    main.run()
