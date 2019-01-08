import argparse
import datetime
import logging
from typing import Any, Dict, List, Optional

import os
import shutil

from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.components.workflows.genedetectionwrapper import GeneDetectionWrapper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources import CSS_STYLE


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
        argument_parser.add_argument('--fasta', help="Input FASTA file", type=str)
        argument_parser.add_argument('--fasta-name', help="Input FASTA file name", type=str)
        argument_parser.add_argument('--fastq-pe', help="Input PE FASTQ files", nargs=2)
        argument_parser.add_argument('--fastq-pe-names', help="Input PE FASTQ file names", nargs=2)
        argument_parser.add_argument('--kmers', help="Kmers to use for assembly")
        argument_parser.add_argument('--trim-reads', help="Perform read trimming", action='store_true')
        argument_parser.add_argument('--database-dir', type=str, required=True)
        argument_parser.add_argument('--detection-method', type=str, choices=['blast', 'srst2'], default='blast')
        argument_parser.add_argument('--report-include-fastq', action='store_true')

        # BLAST specific parameters
        argument_parser.add_argument('--blast-min-percent-identity', type=int, default=90)
        argument_parser.add_argument('--blast-min-percent-coverage', type=int, default=60)

        # SRST2 specific parameters
        argument_parser.add_argument('--srst2-min-cov', type=int, default=90)
        argument_parser.add_argument('--srst2-max-div', type=int, default=10)
        return argument_parser.parse_args()

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        self.__init_report()
        self.__add_analysis_info_section()
        db_data = self.__get_db_metadata()
        if self._args.detection_method == 'blast':
            fasta_file = self.__get_blast_input()
            output = self.__run_gene_detection_blast(fasta_file, db_data)
        else:
            fastq_files = self.__get_srst2_input()
            output = self.__run_gene_detection_srst2(fastq_files, db_data)
        self.__export_output(output)

    def __init_report(self) -> None:
        """
        Initializes the HTML report
        :return: None
        """
        self._report = HtmlReport(self._args.output_html, self._args.output_dir)
        if not os.path.isdir(self._args.output_dir):
            os.makedirs(self._args.output_dir)
        self._report.initialize('Gene detection report', CSS_STYLE)
        self._report.add_pipeline_header(f'Gene detection ({self._args.detection_method})')
        self._report.save()

    def __add_analysis_info_section(self) -> None:
        """
        Adds the report section with the analysis info
        :return: None
        """
        section = HtmlReportSection('Analysis info')
        common_data = [
            ['Analysis date:', datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT)],
            ['Input file(s):', self._helper.determine_input_files(self._args)]
        ]
        if self._args.detection_method == 'blast':
            common_data.extend([
                ['% identity threshold:', f'{self._args.blast_min_percent_identity}%'],
                ['% query covered threshold:', f'{self._args.blast_min_percent_coverage}%']
            ])
        elif self._args.detection_method == 'srst2':
            common_data.extend([
                ['Min. % coverage threshold:', f'{self._args.srst2_min_cov}%'],
                ['Max. % divergence threshold:', f'{self._args.srst2_max_div}%']
            ])
        section.add_table(common_data, table_attributes=[('class', 'information')])
        self._report.add_html_object(section)
        self._report.save()

    def __get_blast_input(self) -> ToolIOFile:
        """
        Returns the input FASTA file
        :return: FASTA file
        """
        if self._args.fasta is not None:
            return ToolIOFile(self._args.fasta)
        else:
            if self._args.trim_reads:
                assembly_input = self._helper.trim_reads(
                    self._args.fastq_pe, self._report, self._args.threads, self._args.report_include_fastq)
            else:
                assembly_input = self._helper.symlink_fastq_pe_input(
                    self._args.fastq_pe, self._args.fastq_pe_names, self._args.working_dir)
            return self._helper.assemble_fastq_reads(assembly_input, self._report, self._args.kmers, self._args.threads)

    def __get_srst2_input(self) -> List[ToolIOFile]:
        """
        Returns the SRST2 input files.
        :return: Paired end FASTQ files
        """
        if self._args.trim_reads is False:
            links = SnakePipelineUtils.symlink_input_files(
                os.path.join(self._args.working_dir, 'input'), self._args.fastq_pe,
                [f'{self._sample_name}_{i}.fastq' for i in (1, 2)])
            return [ToolIOFile(l) for l in links]
        else:
            trimming_output = self._helper.trim_reads(
                self._args.fastq_pe, self._report, self._args.threads, self._args.report_include_fastq)
            return trimming_output.pe

    def __get_db_metadata(self) -> Dict[str, Any]:
        """
        Returns the database information dictionary.
        :return: Database information dictionary
        """
        metadata = {'path': self._args.database_dir}
        if self._args.detection_method == 'blast':
            metadata.update({
                'min_percent_identity': self._args.blast_min_percent_identity,
                'min_coverage': self._args.blast_min_percent_coverage
            })
        return metadata

    def __run_gene_detection_blast(self, fasta_file: ToolIOFile, db_data: Dict[str, Any]) -> \
            GeneDetectionWrapper.GeneDetectionOutput:
        """
        Runs the gene detection workflow.
        :param fasta_file: FASTA file
        :param db_data: Database information dictionary
        :return: None
        """
        wrapper = GeneDetectionWrapper(os.path.join(self._args.working_dir, 'resfinder'))
        wrapper.run_workflow_blast(fasta_file.path, self._sample_name, db_data)
        return wrapper.output

    def __run_gene_detection_srst2(self, fastq_pe: List[ToolIOFile], db_data: Dict[str, Any]) -> \
            GeneDetectionWrapper.GeneDetectionOutput:
        """
        Runs the gene detection workflow in srst2 mode.
        :param fastq_pe: Paired end FASTQ input
        :param db_data: Database information dictionary
        :return: None
        """
        wrapper = GeneDetectionWrapper(os.path.join(self._args.working_dir, 'resfinder'))
        wrapper.run_workflow_srst2([f.path for f in fastq_pe], self._sample_name, db_data)
        return wrapper.output

    def __export_output(self, output: GeneDetectionWrapper.GeneDetectionOutput) -> None:
        """
        Exports the output of the workflow.
        :param output: Output
        :return: None
        """
        self._report.add_html_object(output.report_section)
        output.report_section.copy_files(self._report.output_dir)
        if output.log_file is not None:
            shutil.copyfile(output.log_file, os.path.join(self._report.output_dir, 'log.txt'))
        self._report.save()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main = MainGeneDetection()
    main.run()
