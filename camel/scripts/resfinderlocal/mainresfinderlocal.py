#!/usr/bin/env python
import argparse
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.workflows.genedetectionwrapper import GeneDetectionWrapper, GeneDetectionOutput
from camel.app.components.workflows.readtype import helper_by_read_type
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils


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
        self._sample_name = mainscriptutils.determine_sample_name(self._args)
        self._helper = helper_by_read_type[self._args.read_type](self._args.working_dir, self._sample_name)

    @staticmethod
    def parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        mainscriptutils.add_common_arguments(argument_parser)
        mainscriptutils.add_input_files_arguments(argument_parser)
        mainscriptutils.add_assembly_arguments(argument_parser)
        argument_parser.add_argument('--min-percent-identity', type=int, default=90)
        argument_parser.add_argument('--min-percent-coverage', type=int, default=60)
        argument_parser.add_argument('--resfinder-db', type=str, required=True)
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        # Initialize report
        report = mainscriptutils.init_report(
            self._args.output_html, self._args.output_dir, 'ResFinder local report', 'ResFinder local')
        report.add_html_object(mainscriptutils.generate_analysis_info_section(self._args))
        report.save()

        # Run tools
        fasta_file = self._helper.prepare_fasta_input(report, self._args)
        db_data = self.__get_db_data()
        wrapper = self.__run_gene_detection(fasta_file, db_data, report)

        # Save report
        all_informs = self._helper.informs + [wrapper.informs]
        report.add_html_object(SnakePipelineUtils.create_commands_section(all_informs, self._args.working_dir))
        report.add_html_object(SnakePipelineUtils.create_citations_section(['Zankari_2012-resfinder']))
        report.save()

    def __get_db_data(self) -> Dict[str, Any]:
        """
        Returns the database information dictionary.
        :return: Database information dictionary
        """
        return {
            'path': self._args.resfinder_db,
            'params': {
                'blastn': {
                    'min_percent_identity': self._args.min_percent_identity,
                    'min_coverage': self._args.min_percent_coverage}},
            'metadata': {'name': 'Antibiotic(s)', 'key': 'antibiotics'}
        }

    def __run_gene_detection(self, fasta_file: Path, db_data: Dict[str, Any], report: HtmlReport) -> \
            GeneDetectionOutput:
        """
        Runs the gene detection workflow.
        :param fasta_file: FASTA file
        :param db_data: Database information dictionary
        :return: None
        """
        wrapper = GeneDetectionWrapper(self._args.working_dir / 'resfinder')
        wrapper.run_workflow_blast(fasta_file, self._sample_name, db_data)
        report.add_html_object(wrapper.output.report_section)
        wrapper.output.report_section.copy_files(report.output_dir)
        report.save()
        return wrapper.output


if __name__ == '__main__':
    Camel.get_instance()
    main = MainResFinderLocal()
    main.run()
