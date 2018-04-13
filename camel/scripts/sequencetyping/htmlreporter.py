import os
import shutil

import datetime

from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.resources import CSS_STYLE


class HtmlReporter(object):
    """
    HTML reporter for the sequence typing tool.
    """

    def __init__(self, output_dir, report):
        """
        Initializes the HTML reporter.
        :param output_dir: Output directory
        :param report: Report
        """
        self._output_dir = output_dir
        self._report = report
        self.__create_output_dir()

    def initialize(self, fasta_name, fastq_names):
        """
        Initializes the HTML report.
        :param fasta_name: Name of the FASTA input file (None if FASTQ input)
        :param fastq_names: Name of the FASTQ input files (None if FASTA input)
        :return: None
        """
        # TODO: Fix report
        report = HtmlReport(self._report)
        report.initialize('MLST Output', CSS_STYLE)
        report.add_table([['Analysis date:', datetime.datetime.now().strftime('%d/%m/%Y - %X')],
                          ['Input file(s):', fasta_name if fasta_name is not None else ', '.join(fastq_names)]],
                         table_attributes=[('class', 'information')])
        report.save()

    def __create_output_dir(self):
        """
        Creates the output directory.
        :return: None
        """
        try:
            os.mkdir(self._output_dir)
        except OSError:
            pass

    def add_trimming_section(self, trimming_info):
        """
        Adds the read trimming section to the report.
        :param trimming_info: Trimming info
        :return: None
        """
        report = HtmlReportSection('Read Trimming')
        report.add_header('Read Trimming', 3)
        table_data = [
            ['Input Reads Pairs:', trimming_info['paired_reads_in']],
            ['Both Surviving:', trimming_info['paired_reads_out']],
            ['Forward Only Surviving:', trimming_info['forward_only_reads']],
            ['Reverse Only Surviving:', trimming_info['reverse_only_reads']],
            ['Dropped:', trimming_info['reads_drop']]
        ]
        report.add_table(table_data, table_attributes=[('class', 'information')])

    def add_fastqc_section(self, reports):
        """
        Adds the FastQC section.
        :param reports: FastQC reports
        :return: None
        """
        report = HtmlReportSection("FastQC")
        report.add_header('Post-trimming QC reports', 3)
        output_folder = os.path.join(self._output_dir, 'fastqc')
        os.mkdir(output_folder)
        for fqc_report, orientation in zip(reports, ('forward', 'reverse',)):
            new_name = 'fastqc_{}.html'.format(orientation)
            shutil.copy(fqc_report.path, os.path.join(output_folder, new_name))
            report.add_link_to_file('{}_{}.html'.format(os.path.splitext(fqc_report.basename)[0], orientation),
                                    os.path.join('fastqc', new_name))

    def add_assembly_section(self, sample_name, assembly):
        """
        Adds the assembly section to the report.
        :param sample_name: Sample name
        :param assembly: Assembly FASTA file
        :return: None
        """
        report = HtmlReportSection("Assembly")
        report.add_header('Assembly', 3)
        output_dir = os.path.join(self._output_dir, 'assembly')
        os.mkdir(output_dir)
        assembly_name = 'contigs_{}.fasta'.format(sample_name)
        shutil.copy(assembly.path, os.path.join(output_dir, assembly_name))
        report.add_link_to_file('Assembly (FASTA)', os.path.join('assembly', assembly_name))
