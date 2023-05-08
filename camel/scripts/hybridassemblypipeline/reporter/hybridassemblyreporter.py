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

        self.__add_assemblies_to_download()
        self.report.add_html_object(self._input_informs['trimming_illumina'])
        self.report.add_html_object(self._input_informs['trimming_ont'])
        self.__add_quast_table()
        self.__add_vc_table()
        self.__add_sniffles_table()
        self.__add_mapping_table()

        self.report.add_html_object(self._input_informs['commands'][0].value)
        self.report.save()

    def __add_input_section(self):
        """
        Adds the information about the input data.
        """
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

    def __add_assemblies_to_download(self):
        """
        Adds a new section which allows to download the generated assemblies at different stages of the pipeline.
        """
        fasta_level = ['Flye', 'Medaka', 'POLCA', 'Polypolish', 'Unicycler']
        section = HtmlReportSection('Downloadable assemblies')
        for FASTA in fasta_level:
            fasta_file = Path(self._output_dir) / 'qc' / f'{FASTA}' / 'consensus.fasta'
            relative_path = Path(f'qc/{FASTA}', fasta_file.name)
            section.add_file(fasta_file, relative_path)
            section.add_link_to_file(f'Download {FASTA} assembly (FASTA)', relative_path)
        self.report.add_html_object(section)

    def __add_report_header(self):
        """
        Adds the report header.
        """
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
        """
        Adds the summary QUAST table to the report.
        """
        section = HtmlReportSection('Quast statistics', subtitle='v4.4')
        div_sect = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div = HtmlExpandableDiv('quast_statistics', 'quast')
        div.add_table(self._input_informs['quast'],
                      column_names=['Assembly step', 'N50', 'No of contigs', 'Total length'],
                      table_attributes=[('class', 'data')])
        div_sect.add_html_object(div)
        section.add_html_object(div_sect)
        self.report.add_html_object(section)

    def __add_vc_table(self):
        """
        Adds the variant calling section to the report.
        """
        section = HtmlReportSection('Variant calling statistics', subtitle='Freebayes v1.3.6, Clair3 v1.0.0')
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
        """
        Adds the Sniffles table to the report.
        """
        section = HtmlReportSection('Sniffles statistics', subtitle='v2.0.7')
        div_sect = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div = HtmlExpandableDiv('sniffles_statistics', 'sniffles')
        div.add_table(self._input_informs['sniffles'],
                      column_names=['Assembly step', 'Number of variants', 'Number of indels', 'Number of SNPs'],
                      table_attributes=[('class', 'data')])
        div_sect.add_html_object(div)
        section.add_html_object(div_sect)
        self.report.add_html_object(section)

    def __add_mapping_table(self):
        """
        Adds the mapping statistics to the report.
        """
        section = HtmlReportSection('Mapping statistics', subtitle='BWA v0.7.17, Minimap2 v2.17')
        div_sect = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div = HtmlExpandableDiv('mapping_statistics', 'mapping')
        div.add_table(self._input_informs['mapping'],
                      column_names=['Assembly step', 'Mapping rate (short reads)', 'Median depth (short reads)',
                                    'Mapping rate (long reads)', 'Median depth (long reads)'],
                      table_attributes=[('class', 'data')])
        div_sect.add_html_object(div)
        section.add_html_object(div_sect)
        self.report.add_html_object(section)
