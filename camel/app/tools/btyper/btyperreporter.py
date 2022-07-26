from pathlib import Path
from typing import List, Tuple, Union

import re

from camel.app.camel import Camel
from camel.app.components.html.htmlexpandablediv import HtmlExpandableDiv
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool

class BTyperReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the ResFinder output.
    """

    TITLE = 'BTyper'
    URL_PUBMED = 'https://www.ncbi.nlm.nih.gov/pubmed/{id}'

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the tool.
        :param camel: CAMEL instance
        """
        super().__init__('BTyper Reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError('BTyper input (TSV) is required.')
        if 'btyper' not in self._input_informs:
            raise InvalidInputSpecificationError('BTyper informs are required.')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool
        :return: None
        """
        section = HtmlReportSection(BTyperReporter.TITLE, subtitle=self._input_informs['btyper']['_name'])
        header, data = self.__parse_input_file()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]
        # self._informs['genes'] = data

    def __parse_input_file(self) -> Tuple[List[str], List[List[str]]]:
        """
        Parses the input file.
        :return: Input file header, input file data
        """
        with open(self._tool_inputs['TSV'][0].path) as handle:
            header = handle.readline().strip().split('\t')
            return header, [line.strip().split('\t') for line in handle.readlines()]

    def __generate_output_filename(self) -> str:
        """
        Generates the filename of the tabular output.
        :return: Output filename
        """
        if 'sample_name' in self._parameters:
            return f"btyper-{self._parameters['sample_name'].value}.tsv"
        else:
            return 'btyper.tsv'

    # To Complete "add_output_table"
    # To add: classification, virulence genes, subspecies, taxon name, panC group, clonal complex
    def __add_output_table(
            self, section: HtmlReportSection, header: List[str], data: List[List[Union[str, HtmlTableCell]]]) -> None:
        """
        Adds the output table.
        :param section: Report section
        :param header: Output table header
        :param data: Output table data
        :return: None
        """
        if len(data) > 0:
            div = HtmlExpandableDiv('resfinder_genes', 'genes')
            div.add_table(data, header, [('class', 'data')])
            section.add_html_object(div)
            relative_path = Path('resfinder', self.__generate_output_filename())
            section.add_file(self._tool_inputs['TSV'][0].path, relative_path)
            section.add_link_to_file('Download (TSV)', relative_path)
        else:
            section.add_paragraph('No genes found.')
