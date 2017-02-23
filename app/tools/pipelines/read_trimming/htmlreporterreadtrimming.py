import os

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.tools.export.htmlreporter import HtmlReporter


class HtmlReporterReadTrimming(HtmlReporter):
    """
    Tool to create HTML reports for the read trimming.
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: CAMEL instance
        :return: None
        """
        super(HtmlReporterReadTrimming, self).__init__(camel)
        self.__subfolder = 'quality_control'

    def _create_report(self):
        """
        Creates the HTML report.
        :return: None
        """
        self._report.add_header('Quality Control', 2)
        self.__add_section_qc_pre()
        self.__add_trimming_section()
        self.__add_qc_post_section()
        self._report.add_horizontal_line()

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if 'HTML_Pre' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No pre-trimming reports found")
        if 'FASTQ_PE' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No trimmed paired end reads found")
        if 'HTML_Post' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No post-trimming reports found")
        if 'trimming' not in self._input_informs:
            raise InvalidInputSpecificationError("No trimming info found")
        super(HtmlReporterReadTrimming, self)._check_input()

    def __add_section_qc_pre(self):
        """
        Adds the pre trimming QC section
        :return: None
        """
        self._report.add_header('Pre-trimming', 3)
        for fastqc_report, orientation in zip(self._tool_inputs['HTML_Pre'], ('forward', 'reverse')):
            relative_path = os.path.join(self.__subfolder, 'pre_trimming', 'fastqc_report_{}.html'.format(orientation))
            self._save_file(fastqc_report.path, relative_path)
            self._report.add_link_to_file('FastQC report ({})'.format(orientation), relative_path)

    def __add_qc_post_section(self):
        """
        Adds the post trimming QC section
        :return: None
        """
        self._report.add_header('Post-trimming', 3)
        for fastqc_report, orientation in zip(self._tool_inputs['HTML_Post'], ('forward', 'reverse')):
            relative_path = os.path.join(self.__subfolder, 'post_trimming', 'fastqc_report_{}.html'.format(orientation))
            self._save_file(fastqc_report.path, relative_path)
            self._report.add_link_to_file('FastQC report ({})'.format(orientation), relative_path)

    def __add_trimming_section(self):
        """
        Adds the read trimming section.
        :return: None
        """
        self._report.add_header('Read Trimming', 3)
        trimming_info = self._input_informs['trimming']
        table_data = [
            ['Input Reads Pairs:', trimming_info['paired_reads_in']],
            ['Both Surviving:', trimming_info['paired_reads_out']],
            ['Forward Only Surviving:', trimming_info['forward_only_reads']],
            ['Reverse Only Surviving:', trimming_info['reverse_only_reads']],
            ['Dropped:', trimming_info['reads_drop']]
        ]
        self._report.add_table(table_data, table_attributes=[('class', 'information')])
        self.__add_trimmed_reads()

    def __add_trimmed_reads(self):
        """
        Saves the trimmed reads in the output directory and adds them to the report.
        :return: None
        """
        for trimmed_reads_paired, orientation in zip(self._tool_inputs['FASTQ_PE'], ('forward', 'reverse')):
            filename = 'trimmed_reads_paired_{}.fastq'.format(orientation)
            relative_path = os.path.join(self.__subfolder, 'read_trimming', filename)
            self._save_file(trimmed_reads_paired.path, relative_path)
            self._report.add_link_to_file('Trimmed reads paired ({})'.format(orientation), relative_path)
        self._report.add_line_break()

        if 'FASTQ_SE_FORWARD' in self._tool_inputs:
            relative_path = os.path.join(self.__subfolder, 'read_trimming', 'trimmed_reads_unpaired_forward.fastq')
            self._save_file(self._tool_inputs['FASTQ_SE_FORWARD'][0].path, relative_path)
            self._report.add_link_to_file('Trimmed reads unpaired (forward)', relative_path)
        else:
            self._report.add_link_to_file('Trimmed reads unpaired (forward)', None)

        if 'FASTQ_SE_REVERSE' in self._tool_inputs:
            relative_path = os.path.join(self.__subfolder, 'read_trimming', 'trimmed_reads_unpaired_reverse.fastq')
            self._save_file(self._tool_inputs['FASTQ_SE_REVERSE'][0].path, relative_path)
            self._report.add_link_to_file('Trimmed reads unpaired (reverse)', relative_path)
        else:
            self._report.add_link_to_file('Trimmed reads unpaired (reverse)', None)
