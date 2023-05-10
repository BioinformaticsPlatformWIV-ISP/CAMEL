import datetime
from pathlib import Path
from typing import Any, Dict

import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.tools.tool import Tool
from camel.resources import CSS_STYLE
from camel.resources.javascript import JQUERY_SRC


class HybridAssemblyReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the hybrid assembly pipeline output.
    """

    TITLE = 'Hybrid assembly pipeline'
    MATCH_COLORS = {0: None, 1: 'grey', 2: 'lightgreen', 3: 'green'}
    REPORT_STRUCTURE = [
        ['Read trimming and basic QC', 'trim'],
        ['Quast analysis', 'quast'],
        ['Variant calling analysis', 'vc'],
        ['Sniffles analysis', 'sniffles'],
        ['Mapping analysis', 'mapping']
    ]

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the tool.
        :param camel: CAMEL instance
        """
        super().__init__('Hybrid assembly pipeline reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid.
        :return: None
        """
        if 'TSV_quast' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Combined QUAST input is required ('TSV_quast')")
        if 'quast' not in self._input_informs:
            raise InvalidInputSpecificationError("QUAST informs are required ('quast')")
        if 'TSV_vc' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Combined variant calling input is required ('TSV_vc')")
        if 'freebayes' not in self._input_informs:
            raise InvalidInputSpecificationError("freebayes informs are required ('freebayes')")
        if 'clair3' not in self._input_informs:
            raise InvalidInputSpecificationError("Clair3 informs are required ('clair3')")
        if 'TSV_mapping' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Mapping TSV is required ('TSV_mapping')")
        if 'mapping' not in self._input_informs:
            raise InvalidInputSpecificationError("Mapping informs are required ('mapping')")
        if 'TSV_sniffles' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Sniffles TSV is required ('TSV_sniffles')")
        if 'sniffles' not in self._input_informs:
            raise InvalidInputSpecificationError("Sniffles informs are required ('sniffles')")
        if 'sample_name' not in self._input_informs:
            raise InvalidInputSpecificationError("Sample name is required ('sample_name')")
        if 'pipeline' not in self._input_informs:
            raise InvalidInputSpecificationError("Pipeline informs are required ('pipeline')")
        if 'input' not in self._input_informs:
            raise InvalidInputSpecificationError("Input samples are required ('input')")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._output_dir = Path(self._parameters['output_dir'].value)
        output_html = self._output_dir / 'output.html'
        # Initialize report
        self.report = HtmlReport(output_html, self._output_dir, [JQUERY_SRC])
        self.report.initialize('Hybrid assembly pipeline', CSS_STYLE)
        self.report.add_pipeline_header(HybridAssemblyReporter.TITLE)

        # Input section
        self.__add_input_section()
        self.__add_overview_links()
        self.report.add_module_header('Overview')
        self.__add_overview_section()
        self.report.add_module_header('Read trimming')
        self.report.add_html_object(self._input_informs['trimming_illumina'])
        self.report.add_html_object(self._input_informs['trimming_ont'])
        self.report.add_module_header('Assembly statistics')
        self.__add_quast_section(self._tool_inputs['TSV_quast'][0].path, self._input_informs['quast'])
        self.report.add_module_header('Variant calling')
        self.__add_vc_table(self._tool_inputs['TSV_vc'][0].path)
        self.__add_sniffles_table(self._tool_inputs['TSV_sniffles'][0].path)
        self.__add_mapping_table(self._tool_inputs['TSV_mapping'][0].path)

        self.report.add_html_object(self._input_informs['commands'][0].value)
        self.report.save()

    def __add_input_section(self) -> None:
        """
        Adds the information about the input data.
        :return: None
        """
        input_files = ', '.join([fastq.name for fastq in self._input_informs['input']['illumina']] +
                                [self._input_informs['input']['ont'].name])
        table_data = [
            ['Sample:', self._input_informs['sample_name']],
            ['Analysis date:', datetime.datetime.now()],
            ['Pipeline version:', self._input_informs['pipeline']['version']],
            ['Input files:', input_files]
        ]
        section = HtmlReportSection('Input')
        section.add_table(table_data, table_attributes=[('class', 'information')])
        self.report.add_html_object(section)

    def __add_overview_section(self) -> None:
        """
        Adds a new section which allows to download the generated assemblies at different stages of the pipeline.
        :return: None
        """
        fasta_level = ['Flye', 'Medaka', 'POLCA', 'Polypolish', 'Unicycler']
        section = HtmlReportSection('Overview section')
        assemblies = []
        for fasta_key in fasta_level:
            # TO CHECK ABSOLUTE PATH IN REPORT?
            fasta_file = Path(self._output_dir) / 'qc' / f'{fasta_key}' / 'consensus.fasta'
            relative_path = Path(f'qc/{fasta_key}', fasta_file.name)
            assemblies.append({
                'Assembly step': fasta_key,
                'Download': HtmlTableCell('Download (FASTA)', link=str(relative_path))})
        df = pd.DataFrame(assemblies)
        # noinspection PyTypeChecker
        section.add_table(
            list(df.itertuples(index=False, name=None)),
            column_names=df.columns,
            table_attributes=[('class', 'data')]
        )
        section.add_paragraph('The long-read-first assembly consists of four main steps: \n'
                              '1) Quality control and pre-processing of the long and short reads\n'
                              '2) Assembling the long reads using Flye and polishing with Medaka\n'
                              '3) Polishing using short reads with POLCA, then Polypolish\n'
                              '4) Quality assessment using QUAST and variant callers\n\n'
                              'The short-read-first assembly consists of running Unicycler using the same input'
                              ' as the long-read-first approach.')
        self.report.add_html_object(section)

    def __add_overview_links(self) -> None:
        """
        Adds the report header.
        :return: None
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

    def __add_quast_section(self, path_tsv: Path, quast_informs: Dict[str, Any]) -> None:
        """
        Adds the summary QUAST table to the report.
        :param path_tsv: Path to the combined QUAST TSV file
        :param quast_informs: QUAST informs
        :return: None
        """
        section = HtmlReportSection('QUAST statistics', subtitle=quast_informs['_name'])
        data_quast = pd.read_table(path_tsv)
        # noinspection PyTypeChecker
        section.add_table(
            list(data_quast.itertuples(index=False, name=None)),
            column_names=data_quast.columns,
            table_attributes=[('class', 'data')])
        self.report.add_html_object(section)

    def __add_vc_table(self, path_tsv: Path) -> None:
        """
        Adds the variant calling section to the report.
        :path: Path to the tsv containing variant calling statistics
        :return: None
        """
        subtitle = ', '.join([self._input_informs[x]['_name'] for x in ('freebayes', 'clair3')])
        section = HtmlReportSection('Variant calling (short reads)', subtitle=subtitle)
        data_vc = pd.read_table(path_tsv)
        # noinspection PyTypeChecker
        section.add_table(
            list(data_vc.itertuples(index=False, name=None)),
            column_names=data_vc.columns,
            table_attributes=[('class', 'data')])
        self.report.add_html_object(section)

    def __add_sniffles_table(self, path_tsv: Path) -> None:
        """
        Adds the Sniffles table to the report.
        :path: Path to the tsv containing variant calling from sniffles
        :return: None
        """
        section = HtmlReportSection('Sniffles statistics', subtitle=self._input_informs['sniffles']['_name'])
        data_sniffles = pd.read_table(path_tsv)
        section.add_table(
            list(data_sniffles.itertuples(index=False, name=None)),
            column_names=data_sniffles.columns,
            table_attributes=[('class', 'data')])
        self.report.add_html_object(section)

    def __add_mapping_table(self, path_tsv: Path) -> None:
        """
        Adds the mapping statistics to the report.
        :path_tsv: Path to the tsv containing mapping statistics
        :return: None
        """
        section = HtmlReportSection('Mapping statistics', subtitle=self._input_informs['mapping']['_name'])
        data_mapping = pd.read_table(path_tsv)
        section.add_table(
            list(data_mapping.itertuples(index=False, name=None)),
            column_names=data_mapping.columns,
            table_attributes=[('class', 'data')]
        )
        self.report.add_html_object(section)
