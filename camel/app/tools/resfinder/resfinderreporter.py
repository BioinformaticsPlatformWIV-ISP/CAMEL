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

class ResFinderReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the ResFinder output.
    """

    TITLE = 'ResFinder'
    URL_PUBMED = 'https://www.ncbi.nlm.nih.gov/pubmed/{id}'

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the tool.
        :param camel: CAMEL instance
        """
        super().__init__('ResFinder Reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError('ResFinder input (TSV) is required.')
        if 'resfinder' not in self._input_informs:
            raise InvalidInputSpecificationError('ResFinder informs are required.')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool
        :return: None
        """
        section = HtmlReportSection(ResFinderReporter.TITLE, subtitle=self._input_informs['resfinder']['_name'])
        header, data = self.__parse_input_file()
        data = self.__add_pubmed_links(data)
        # section.add_paragraph(f"Database: <i>{self._input_informs['resfinder']['database']}</i>")
        # self.__add_output_table(section, header, data)
        # section.add_paragraph('Last update: {}'.format(self._input_informs['resfinder']['last_update']))
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]
        self._informs['genes'] = data

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
            return f"resfinder-{self._parameters['sample_name'].value}.tsv"
        else:
            return 'resfinder.tsv'

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

    def __add_pubmed_links(self, data: List[List[str]]) -> List[List[Union[str, HtmlTableCell]]]:
        """
        Adds PubMed links for mutations that have an associated PMID.
        :param data: Data
        :return: Data with links added
        """
        print(data)
        edited_data = []
        n_fields = len(data[0])
        for i in range(1, len(data)):
            row = data[i]
            if row[n_fields-1]:
                pmid = row[n_fields-1]
                row[n_fields-1] = HtmlTableCell(str(pmid), link=ResFinderReporter.URL_PUBMED.format(id=pmid))
            edited_data.append(row)
        return edited_data

