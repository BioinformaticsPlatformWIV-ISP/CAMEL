#!/usr/bin/env python
import argparse
import logging
from pathlib import Path
from typing import Optional, List, Dict

import yaml

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.components.workflows.genedetectionwrapper import GeneDetectionWrapper
from camel.app.io.tooliofile import ToolIOFile


class MainSCCmecFinder(object):
    """
    This tool is used to run the SCCmecFinder tool.
    """

    def __init__(self, args: Optional[argparse.Namespace] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments, if not set they are removed from the command line
        """
        self._camel = Camel()
        self._args = MainSCCmecFinder._parse_arguments() if args is None else args
        self._sample_name = MainScriptHelper.determine_sample_name(self._args)
        self._helper = MainScriptHelper(self._args.working_dir, self._sample_name)
        self._report = None

    @staticmethod
    def _parse_arguments() -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        MainScriptHelper.add_common_arguments(argument_parser)
        MainScriptHelper.add_assembly_arguments(argument_parser)
        MainScriptHelper.add_input_files_arguments(argument_parser)
        argument_parser.add_argument('--db-mec-genes', help="Database containing mec genes.", required=True)
        argument_parser.add_argument('--profiles-mec-genes', help="Profiles for the mec genes", required=True)
        return argument_parser.parse_args()

    @staticmethod
    def __get_matching_complex(detected_genes: List[str], genes_by_complex: Dict[str, List[str]]) -> \
            Optional[str]:
        """
        Returns the matching complex (if there is one).
        :param genes_by_complex: Genes by complex
        :return: Complex (or None if there is none found)
        """
        for complex_, genes in genes_by_complex.items():
            if all(g in detected_genes for g in genes):
                return complex_
        logging.debug("No complex found")

    def run(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._report = self._helper.init_report(
            self._args.output_html, self._args.output_dir, 'SCCmecFinder ouptput', 'SCCmecFinder (local)')
        self._helper.export_analysis_info_section(self._report, self._helper.determine_input_files(self._args))
        input_files = self._helper.symlink_input_files(self._args.fasta, self._args.fastq_pe)
        fasta_file = self._helper.get_blast_input(input_files, self._args, self._report)
        detected_genes = self.__run_blast(fasta_file)
        self.__add_mec_type_overview(detected_genes)

    def __run_blast(self, fasta_file: ToolIOFile) -> List[str]:
        """
        Runs BLAST on the mec genes database.
        :param fasta_file: Input FASTA file
        :return: List of detected genes
        """
        wrapper = GeneDetectionWrapper(str(Path(self._args.working_dir) / 'meca'))
        wrapper.run_workflow_blast(
          fasta_file.path, self._sample_name, {'path': self._args.db_mec_genes}, self._args.threads)
        self._report.add_html_object(wrapper.output.report_section)
        wrapper.output.report_section.copy_files(self._report.output_dir)
        self._report.save()
        return [d.locus.split(':')[0] for d in wrapper.output.detected_hits]

    def __add_mec_type_overview(self, detected_genes: List[str]) -> None:
        """
        Determines the mec type based on the detected genes and adds it to the report.
        :param detected_genes: Detected genes
        :return: None
        """
        with open(self._args.profiles_mec_genes) as handle:
            profiles = yaml.load(handle)
        ccr_complex = MainSCCmecFinder.__get_matching_complex(detected_genes, profiles['ccr_genes_complexes'])
        mec_complex = MainSCCmecFinder.__get_matching_complex(detected_genes, profiles['mec_genes_complexes'])
        sccmec_type = MainSCCmecFinder.__get_matching_complex(detected_genes, profiles['SCC_mec_types'])
        section = HtmlReportSection('SCCmec type', 3)
        section.add_table([
            ['SCCmec type:', sccmec_type if sccmec_type is not None else '-'],
            ['<i>mec</i> class:', mec_complex if mec_complex is not None else '-'],
            ['<i>ccr</i> class:', ccr_complex if ccr_complex is not None else '-']
        ], None, [('class', 'information')])
        self._report.add_html_object(section)
        self._report.save()


if __name__ == '__main__':
    c = Camel.get_instance()
    main = MainSCCmecFinder()
    main.run()
