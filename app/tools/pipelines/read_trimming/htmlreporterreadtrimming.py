import os

from app.components.html.htmlreportsection import HtmlReportSection
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliovalue import ToolIOValue
from app.tools.tool import Tool


class HtmlReporterReadTrimming(Tool):
    """
    Tool to create HTML reports for the read trimming.
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('HTML Reporter', '0.1', camel)
        self.__subfolder = 'read_trimming'
        self._report_section = None

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self._report_section = HtmlReportSection('Read Trimming')
        self.__add_qc_pre_section()
        self.__add_trimming_section()
        self.__add_qc_post_section()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section, False)]

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
        super()._check_input()

    def __add_qc_pre_section(self):
        """
        Adds the pre trimming QC section
        :return: None
        """
        self._report_section.add_header('Pre-trimming', 3)
        for fastqc_report, orientation in zip(self._tool_inputs['HTML_Pre'], ('forward', 'reverse')):
            relative_path = os.path.join(self.__subfolder, 'pre_trimming', 'fastqc_report_{}.html'.format(orientation))
            self._report_section.add_file(fastqc_report.path, relative_path)
            self._report_section.add_link_to_file('FastQC report ({})'.format(orientation), relative_path)

    def __add_qc_post_section(self):
        """
        Adds the post trimming QC section
        :return: None
        """
        self._report_section.add_header('Post-trimming', 3)
        for fastqc_report, orientation in zip(self._tool_inputs['HTML_Post'], ('forward', 'reverse')):
            relative_path = os.path.join(self.__subfolder, 'post_trimming', 'fastqc_report_{}.html'.format(orientation))
            self._report_section.add_file(fastqc_report.path, relative_path)
            self._report_section.add_link_to_file('FastQC report ({})'.format(orientation), relative_path)

    def __add_trimming_section(self):
        """
        Adds the read trimming section.
        :return: None
        """
        self._report_section.add_header('Read Trimming', 3)
        trimming_info = self._input_informs['trimming']
        table_data = [
            ['Input Reads Pairs:', trimming_info['paired_reads_in']],
            ['Both Surviving:', trimming_info['paired_reads_out']],
            ['Forward Only Surviving:', trimming_info['forward_only_reads']],
            ['Reverse Only Surviving:', trimming_info['reverse_only_reads']],
            ['Dropped:', trimming_info['reads_drop']]
        ]
        self._report_section.add_table(table_data, table_attributes=[('class', 'information')])
        self.__add_trimmed_reads()

    def __add_trimmed_reads(self):
        """
        Saves the trimmed reads in the output directory and adds them to the report.
        :return: None
        """
        for trimmed_reads_paired, orientation in zip(self._tool_inputs['FASTQ_PE'], ('forward', 'reverse')):
            filename = 'trimmed_reads_paired_{}.fastq'.format(orientation)
            relative_path = os.path.join(self.__subfolder, 'read_trimming', filename)
            self._report_section.add_file(trimmed_reads_paired.path, relative_path)
            self._report_section.add_link_to_file('Trimmed reads paired ({})'.format(orientation), relative_path)
        self._report_section.add_line_break()

        if 'FASTQ_SE_FORWARD' in self._tool_inputs:
            relative_path = os.path.join(self.__subfolder, 'read_trimming', 'trimmed_reads_unpaired_forward.fastq')
            self._report_section.add_file(self._tool_inputs['FASTQ_SE_FORWARD'][0].path, relative_path)
            self._report_section.add_link_to_file('Trimmed reads unpaired (forward)', relative_path)
        else:
            self._report_section.add_link_to_file('Trimmed reads unpaired (forward)', None)

        if 'FASTQ_SE_REVERSE' in self._tool_inputs:
            relative_path = os.path.join(self.__subfolder, 'read_trimming', 'trimmed_reads_unpaired_reverse.fastq')
            self._report_section.add_file(self._tool_inputs['FASTQ_SE_REVERSE'][0].path, relative_path)
            self._report_section.add_link_to_file('Trimmed reads unpaired (reverse)', relative_path)
        else:
            self._report_section.add_link_to_file('Trimmed reads unpaired (reverse)', None)
