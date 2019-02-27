from typing import List, Tuple, Union

import os
import re

from camel.app.camel import Camel
from camel.app.components.html.htmlexpandablediv import HtmlExpandableDiv
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class PointFinderReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the PointFinder output.
    """

    TITLE = 'PointFinder'
    URL_PUBMED = 'https://www.ncbi.nlm.nih.gov/pubmed/{id}'

    def __init__(self, camel: Camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('PointFinder Reporter', '0.1', camel)
        self._section = HtmlReportSection(PointFinderReporter.TITLE)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Tabular PointFinder input (TSV) is required')
        if 'pointfinder' not in self._input_informs:
            raise InvalidInputSpecificationError('PointFinder informs are required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        header, data = self.__parse_input_file()
        data = self.__add_pubmed_links(data)
        self._section.add_paragraph(f"Database: <i>{self._input_informs['pointfinder']['database']}</i>")
        self.__add_output_table(header, data)
        self._section.add_paragraph('Last update: {}'.format(self._input_informs['pointfinder']['last_update']))
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]
        self._informs['mutations'] = data

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
            return f"pointfinder-{self._parameters['sample_name'].value}.tsv"
        else:
            return 'pointfinder.tsv'

    def __add_output_table(self, header: List[str], data: List[List[Union[str, HtmlTableCell]]]) -> None:
        """
        Adds the output table.
        :param header: Output table header
        :param data: Output table data
        :return: None
        """
        div = HtmlExpandableDiv('pointfinder_mutations', 'mutations')
        div.add_table(data, header, [('class', 'data')])
        self._section.add_html_object(div)
        relative_path = os.path.join('pointfinder', self.__generate_output_filename())
        self._section.add_file(self._tool_inputs['TSV'][0].path, relative_path)
        self._section.add_link_to_file('Download (TSV)', relative_path)

    def __add_pubmed_links(self, data: List[List[str]]) -> List[List[Union[str, HtmlTableCell]]]:
        """
        Adds PubMed links for mutations that have an associated PMID.
        :param data: Data
        :return: Data with links added
        """
        edited_data = []
        for row in data:
            m = re.match('^(\\d+)$', row[-1])
            if not m:
                new_row = row
            else:
                pmid = int(m.group(1))
                new_row = row[:-1]
                # noinspection PyTypeChecker
                new_row.append(HtmlTableCell(str(pmid), link=PointFinderReporter.URL_PUBMED.format(id=pmid)))
            edited_data.append(new_row)
        return edited_data
