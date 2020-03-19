#!/usr/bin/env python
import argparse
import datetime
import logging
from typing import Any, Dict, Optional, Sequence

import os

from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.components.workflows.genedetectionwrapper import GeneDetectionWrapper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources import CSS_STYLE


class MainResFinderLocal(object):
    """
    This class is used to run the main ResFinder local script.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments (optional)
        """
        self._args = MainResFinderLocal.parse_arguments(args)
        self._sample_name = MainScriptHelper.determine_sample_name(self._args)
        self._helper = MainScriptHelper(self._args.working_dir, self._sample_name)
        self._report = None

    @staticmethod
    def parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        MainScriptHelper.add_common_arguments(argument_parser)
        argument_parser.add_argument('--fasta', help="Input FASTA file")
        argument_parser.add_argument('--fasta-name', help="Input FASTA file name")
        argument_parser.add_argument('--fastq-pe', help="Input PE FASTQ files", nargs=2)
        argument_parser.add_argument('--fastq-pe-names', help="Input PE FASTQ file names", nargs=2)
        argument_parser.add_argument('--trim-reads', help="Perform read trimming", action='store_true')
        argument_parser.add_argument('--kmers', help="Kmers to use for assembly")
        argument_parser.add_argument('--min-percent-identity', type=int, default=90)
        argument_parser.add_argument('--min-percent-coverage', type=int, default=60)
        argument_parser.add_argument('--resfinder-db', type=str, required=True)
        argument_parser.add_argument('--report-include-fastq', action='store_true')
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        self.__init_report()
        self.__add_analysis_info_section()
        input_files = self._helper.symlink_input_files(self._args.fasta, self._args.fastq_pe)
        fasta_file = self._helper.get_blast_input(input_files, self._args, self._report)
        db_data = self.__get_db_data()
        self.__run_gene_detection(fasta_file, db_data)

    def __init_report(self) -> None:
        """
        Initializes the HTML report
        :return: None
        """
        self._report = HtmlReport(self._args.output_html, self._args.output_dir)
        if not os.path.isdir(self._args.output_dir):
            os.makedirs(self._args.output_dir)
        self._report.initialize('ResFinder local report', CSS_STYLE)
        self._report.add_pipeline_header('ResFinder local')

    def __add_analysis_info_section(self) -> None:
        """
        Adds the report section with the analysis info
        :return: None
        """
        section = HtmlReportSection('Analysis info')
        section.add_table([
            ['Analysis date:', datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT)],
            ['Input file(s):', self._helper.determine_input_files(self._args)],
            ['Selected % identity threshold:', f'{self._args.min_percent_identity}%'],
            ['Selected % query covered threshold:', f'{self._args.min_percent_coverage}%']
        ], table_attributes=[('class', 'information')])
        self._report.add_html_object(section)
        self._report.save()

    def __get_db_data(self) -> Dict[str, Any]:
        """
        Returns the database information dictionary.
        :return: Database information dictionary
        """
        return {
            'path': self._args.resfinder_db,
            'min_percent_identity': self._args.min_percent_identity,
            'min_coverage': self._args.min_percent_coverage,
            'extra_column': {'name': 'Antibiotic(s)', 'key': 'antibiotics'}
        }

    def __run_gene_detection(self, fasta_file: ToolIOFile, db_data: Dict[str, Any]) -> None:
        """
        Runs the gene detection workflow.
        :param fasta_file: FASTA file
        :param db_data: Database information dictionary
        :return: None
        """
        wrapper = GeneDetectionWrapper(os.path.join(self._args.working_dir, 'resfinder'))
        wrapper.run_workflow_blast(fasta_file.path, self._sample_name, db_data)
        self._report.add_html_object(wrapper.output.report_section)
        wrapper.output.report_section.copy_files(self._report.output_dir)
        self._report.save()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main = MainResFinderLocal()
    main.run()
