import datetime
from pathlib import Path

from camel.app.camel import Camel
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlexpandablediv import HtmlExpandableDiv
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.tools.tool import Tool


class HybridAssemblyReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the hybrid assembly pipeline output.
    """

    TITLE = 'Hybrid assembly pipeline'
    MATCH_COLORS = {0: None, 1: 'grey', 2: 'lightgreen', 3: 'green'}
    REPORT_STRUCTURE = [['Read trimming and basic QC', 'trim'], ['Quast analysis', 'quast'],
                        ['Variant calling analysis', 'vc'], ['Sniffles analysis', 'sniffles'],
                        ['Mapping analysis', 'mapping']]

    def __init__(self, camel: Camel, output_directory: Path) -> None:
        """
        Initializes the tool.
        :param camel: CAMEL instance
        """
        self._output_dir = output_directory
        self._output_html = self._output_dir / 'output.html'
        super().__init__('Hybrid assembly pipeline reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid.
        :return: None
        """
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.report = HtmlReport(self._output_html, self._output_dir)

        self._check_input()

        self.report.add_pipeline_header(HybridAssemblyReporter.TITLE)
        self.__add_input_section()

        self.__add_report_header()

        self.report.add_html_object(self._input_informs['trimming_illumina'])
        self.__add_quast_table()
        self.__add_vc_table()
        self.__add_sniffles_table()
        self.__add_mapping_table()

        self.report.add_html_object(self._input_informs['commands'][0].value)
        self.report.save()

    def __add_input_section(self):
        table_data = [
            ['Sample:', self._input_informs['sample_name']],
            ['Analysis date:', datetime.datetime.now()],
            ['Pipeline version:', self._input_informs['pipeline']['version']],
            ['Input files:', [*[i for i in self._input_informs['fastq_input']],
                              self._input_informs['fastq_se_input']]],
        ]
        section = HtmlReportSection('Input')
        section.add_table(table_data, table_attributes=[('class', 'information')])
        self.report.add_html_object(section)

    def __add_report_header(self):
        self.report.add_module_header('Sections')
        section = HtmlReportSection(None)

        overview_list = HtmlElement('ul')
        for title, key in HybridAssemblyReporter.REPORT_STRUCTURE:
            list_item = HtmlElement('li')
            list_item.add_html_object(HtmlElement('a', title, [('href', '#{}'.format(key))]))
            overview_list.add_html_object(list_item)
        section.add_html_object(overview_list)
        self.report.add_html_object(section)

    def __add_quast_table(self):
        section = HtmlReportSection('Quast statistics', subtitle='0.1')
        div_sect = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div = HtmlExpandableDiv('quast_statistics', 'quast')
        div.add_table(self._input_informs['quast'],
                      column_names=['Assembly step', 'N50', 'No of contigs', 'Total length'],
                      table_attributes=[('class', 'data')])
        div_sect.add_html_object(div)
        section.add_html_object(div_sect)
        self.report.add_html_object(section)

    def __add_vc_table(self):
        section = HtmlReportSection('Variant calling statistics', subtitle='0.1')
        div_sect = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div = HtmlExpandableDiv('variant_calling_statistics', 'vc')
        div.add_table(self._input_informs['vc'],
                      column_names=['Assembly step', 'Freebayes total variants' 'Freebayes indels', 'Freebayes SNPs',
                                    'Clair3 total variant', 'Clair3 indels', 'Clair3 SNPs'],
                      table_attributes=[('class', 'data')])
        div_sect.add_html_object(div)
        section.add_html_object(div_sect)
        self.report.add_html_object(section)

    def __add_sniffles_table(self):
        section = HtmlReportSection('Sniffles statistics', subtitle='0.1')
        div_sect = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div = HtmlExpandableDiv('sniffles_statistics', 'sniffles')
        div.add_table(self._input_informs['sniffles'],
                      column_names=['Assembly step', 'Number of variants', 'Number of indels', 'Number of SNPs'],
                      table_attributes=[('class', 'data')])
        div_sect.add_html_object(div)
        section.add_html_object(div_sect)
        self.report.add_html_object(section)

    def __add_mapping_table(self):
        section = HtmlReportSection('Mapping statistics', subtitle='0.1')
        div_sect = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div = HtmlExpandableDiv('mapping_statistics', 'mapping')
        div.add_table(self._input_informs['mapping'],
                      column_names=['Assembly step', 'Mapping rate (short reads)', 'Median depth (short reads)',
                                    'Mapping rate (long reads)', 'Median depth (long reads)'],
                      table_attributes=[('class', 'data')])
        div_sect.add_html_object(div)
        section.add_html_object(div_sect)
        self.report.add_html_object(section)
