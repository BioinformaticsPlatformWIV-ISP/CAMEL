from pathlib import Path
from typing import List, Tuple, Union

from camel.app.camel import Camel
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
    URL_NUCCORE = 'https://www.ncbi.nlm.nih.gov/nuccore/{id}'
    URL_PUBMED = 'https://pubmed.ncbi.nlm.nih.gov/{id}'

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
        data = self.__add_links_and_format_data(data)
        self.__add_output_table(section, header, data)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]
        self._informs['genes'] = data

    def __parse_input_file(self) -> Tuple[List[List[str]], List[List[List[str]]]]:
        """
        Parses the input file.
        :return: Input file header, input file data
        """
        header, data = [], []
        for input_file in self._tool_inputs['TSV']:
            with open(input_file.path) as handle:
                header.append(handle.readline().strip().split('\t'))
                data.append([line.strip().split('\t') for line in handle.readlines()])
        return header, data

    def __generate_output_filename(self, tool_name: str) -> str:
        """
        Generates the filename of the tabular output.
        :return: Output filename
        """
        if 'sample_name' in self._parameters:
            return f"{tool_name}-{self._parameters['sample_name'].value}.tsv"
        else:
            return f'{tool_name}.tsv'

    def __add_output_table(
            self, section: HtmlReportSection, header: List[List[str]],
            data: List[List[Union[str, HtmlTableCell]]]) -> None:
        """
        Adds the output table.
        :param section: Report section
        :param header: Output table header
        :param data: Output table data
        :return: None
        """
        tool_ids = ['resfinder', 'pointfinder']
        for i in range(len(data)):
            table = data[i]
            if len(table) > 0:
                relative_path = Path('resfinder', self.__generate_output_filename(tool_ids[i]))
                section.add_header(f'{tool_ids[i].upper()} Results', 4)
                section.add_table(table, header[i], [('class', 'data')])
                section.add_file(self._tool_inputs['TSV'][i].path, relative_path)
                section.add_link_to_file('Download (TSV)', relative_path)
            else:
                section.add_paragraph('No genes found.')

    def __add_links_and_format_data(self, data: List[List[List[str]]]) -> List[List[Union[str, HtmlTableCell]]]:
        """
        Adds PubMed links for mutations that have an associated PMID.
        :param data: Data
        :return: Data with links added
        """
        table_results = []
        for table in data:
            edited_data = []
            n_fields = len(table[0])
            for i in range(len(table)):
                row = table[i]
                if row[n_fields - 1] is not None:
                    access = row[n_fields - 1]
                    if table_results == []:
                        publi = HtmlTableCell(str(access), link=ResFinderReporter.URL_NUCCORE.format(id=access))
                    else:
                        publi = HtmlTableCell(str(access), link=ResFinderReporter.URL_PUBMED.format(id=access))
                    row[n_fields - 1] = publi
                    try:
                        row[3] = f'{eval(row[3]):.2f}'
                    except NameError:
                        pass
                edited_data.append(row)
            table_results.append(edited_data)
        return table_results
