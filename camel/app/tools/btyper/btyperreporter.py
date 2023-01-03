from pathlib import Path
from typing import List, Tuple, Union
import re

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class BTyperReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the BTyper output.
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
        Executes the BTyper reporter tool
        :return: None
        """
        section = HtmlReportSection(BTyperReporter.TITLE, subtitle=self._input_informs['btyper']['_name'])

        # Add output tables
        headers, data = self.__parse_input_file()
        self.__add_output_table(section, headers[0], data[0], 'Identification (ANI)')
        self.__add_output_table(section, headers[1], data[1], 'Virulence genes')
        self.__add_output_table(section, headers[2], data[2], 'Typing')

        # Add link to TSV file
        relative_path = Path('btyper') / self._tool_inputs['TSV'][0].path.name
        section.add_link_to_file('Download complete results (TSV)', relative_path)
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)

        # Tool output
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]

    def __sanitize_table_line(self, input_list: List[str]) -> List[str]:
        """
        Format nicely the entries in the HTML output table
        :input_list: split line from the raw output table
        :return: List
        """
        import pprint
        pprint.pprint(input_list)
        for k in range(len(input_list)):
            if 'not performed' in input_list[k]:
                input_list[k] = 'Analysis not performed'

            # Remove empty brackets
            input_list[k] = re.sub(r'\(\)', '', input_list[k])

            # Add space before brackets
            input_list[k] = re.sub(r'(\w)\(', r'\1 (', input_list[k])

            #  Italics for gene names
            match_re = re.findall(r'[a-z]{3}[A-Z]', input_list[k])
            if match_re:
                input_list[k] = '{} ({})'.format(
                    input_list[k][:3], '; '.join(['<i>{}</i>'.format(gene) for gene in match_re]))
        return input_list

    def __parse_input_file(self) -> Tuple[List[List[str]], List[List[List[str]]]]:
        """
        Parses the input file.
        :return: Input file header, input file data
        """
        with open(self._tool_inputs['TSV'][0].path) as handle:
            handle.readline()
            header_part1 = ['Species', 'Sub-species']
            header_part2 = ['Anthrax toxin genes', 'Emetic toxin genes', 'Diarrheal toxin <b>Nhe</b> genes',
                            'Diarrheal toxin <b>Hbl</b> genes', 'Diarrheal toxin <b>CytK</b>',
                            'Sphingomyelinase <b>Sph</b>', 'Capsule <b>Cap</b> genes', 'Capsule <b>Has</b> genes',
                            'Capsule <b>Bps</b> genes', '<b>Bt</b> toxin genes']
            header_part3 = ['PubMLST ST', 'Adjusted <i>panC</i> group', 'Final taxon name']
            header_table = [header_part1, header_part2, header_part3]
            output_table = [[], [], []]
            for line in handle.readlines():
                spl = self.__sanitize_table_line(line.strip().split('\t'))
                output_table[0].append(spl[2:4])
                output_table[1].append(spl[6:16])
                output_table[2].append(spl[16:19])
            return header_table, output_table

    def __generate_output_filename(self, prefix: str) -> str:
        """
        Generates the filename of the tabular output.
        :return: Output filename
        """
        if 'sample_name' in self._parameters:
            return f"btyper-{self._parameters['sample_name'].value}_{prefix}.tsv"
        else:
            return f'btyper_{prefix}.tsv'

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
            section.add_header(f'{prefix} analysis', level=3)
            section.add_table(data, header, [('class', 'data')])
        else:
            section.add_paragraph('No results.')
