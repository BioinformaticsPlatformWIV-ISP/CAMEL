from pathlib import Path
from typing import List, Tuple, Union

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class FastANIReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the FastANI output.
    """

    TITLE = 'FastANI'
    URL_PUBMED = 'https://www.ncbi.nlm.nih.gov/pubmed/{id}'

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the tool.
        :param camel: CAMEL instance
        """
        super().__init__('FastANI Reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FastANI input (TSV) is required.')
        if 'fastani' not in self._input_informs:
            raise InvalidInputSpecificationError('FastANI informs are required.')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the FastANI reporter tool
        :return: None
        """
        section = HtmlReportSection(FastANIReporter.TITLE, subtitle=self._input_informs['fastani']['_name'])

        # Add output tables
        header, data = self.__parse_input_file()
        self.__add_output_table(section, header, data)

        # Add link to TSV file
        relative_path = Path('fastani') / self._tool_inputs['TSV'][0].path.name
        section.add_link_to_file('Download complete results (TSV)', relative_path)
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    def __parse_input_file(self) -> Tuple[List[str], List[List[str]]]:
        """
        Parses the FastANI input file.
        :return: Input file header, input file data
        """
        with open(self._tool_inputs['TSV'][0].path) as handle:
            header = ['Query', 'ANI', 'Orthologous matches', 'Total seq fragments']
            output_table = []
            for line in handle.readlines():
                output_table.append(self.__format_output_table_line(line.strip().split()[1:]))
            return header, output_table

    def __format_output_table_line(self, table_line: List[str]) -> List[str]:
        """
        Formats a line in the output table to have only filenames and 2 significant digits on ANI.
        :table_line: input split line from the FastANI output table.
        :return: formatted split line
        """

        # Remove directories from filenames
        table_line[0] = Path(table_line[0]).name

        # Format ANI percentage to two significant digits
        table_line[1] = '{:.2f}'.format(float(table_line[2]))
        return table_line

    def __generate_output_filename(self) -> str:
        """
        Generates the filename of the tabular output.
        :return: Output filename
        """
        if 'sample_name' in self._parameters:
            return f"fastani-{self._parameters['sample_name'].value}.tsv"
        else:
            return f'fastani.tsv'

    def __add_output_table(
            self, section: HtmlReportSection, header: List[str],
            data: List[List[Union[str, HtmlTableCell]]]) -> None:
        """
        Adds an output table to the HTML report.
        :param section: Report section
        :param header: Output table header
        :param data: Output table data
        :return: None
        """
        if len(data) > 0:
            section.add_header(f'FastANI results', level=3)
            section.add_table(data, header, [('class', 'data')])
        else:
            section.add_paragraph('No results.')
