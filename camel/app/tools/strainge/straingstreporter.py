from pathlib import Path
from typing import List, Tuple, Union

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class StrainGSTReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the StrainGST output.
    """

    TITLE = 'StrainGST'
    URL_PUBMED = 'https://www.ncbi.nlm.nih.gov/pubmed/{id}'

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the tool.
        :param camel: CAMEL instance
        """
        super().__init__('StrainGST Reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError('StrainGST strains input (TSV) is required.')
        if 'straingst' not in self._input_informs:
            raise InvalidInputSpecificationError('StrainGST informs are required.')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the StrainGST reporter tool
        :return: None
        """
        section = HtmlReportSection(StrainGSTReporter.TITLE, subtitle=self._input_informs['straingst']['_name'])

        suffix_read_type = self._parameters["suffix"].value.capitalize()
        if suffix_read_type == 'Ont':
            suffix_read_type = 'Nanopore'

        # Add output tables
        header, data = self.__parse_input_file()
        self.__add_output_table(section, header, data,
                                f'StrainGST strain identification - {suffix_read_type}')

        # Add link to TSV file
        relative_path = Path('straingst') / self._tool_inputs['TSV'][0].path.name
        section.add_link_to_file('Download complete results (TSV)', relative_path)
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)

        # Tool output
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]

    def __parse_input_file(self) -> Union[bool, Tuple[List[str], List[List[str]]]]:
        """
        Parses the input file.
        :return: Input file header, input file data
        """
        with open(self._tool_inputs['TSV'][0].path) as handle:
            handle.readline()
            header = ['Strain', 'Coverage', 'Evenness', 'Relative abundance', 'Score']
            output_table = []
            for line in handle.readlines():
                spl = line.strip().split('\t')
                output_table.append([spl[1], spl[5], spl[9], spl[11], spl[14]])
        return header, output_table

    def __generate_output_filename(self, prefix: str) -> str:
        """
        Generates the filename of the tabular output.
        :return: Output filename
        """
        if 'sample_name' in self._parameters:
            return f"straingst-{self._parameters['sample_name'].value}_{prefix}.tsv"
        else:
            return f'straingst_{prefix}.tsv'

    def __add_output_table(
            self, section: HtmlReportSection, header: List[str],
            data: List[List[Union[str, HtmlTableCell]]], prefix: str) -> None:
        """
        Adds an output table to the HTML report.
        :param section: Report section
        :param header: Output table header
        :param data: Output table data
        :param prefix: Prefix for the table to add
        :return: None
        """
        if len(data) > 0:
            section.add_header(f'{prefix}', level=3)
            section.add_table(data, header, [('class', 'data')])
        else:
            section.add_paragraph('No results.')
