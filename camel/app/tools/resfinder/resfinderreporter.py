from pathlib import Path
from typing import List, Tuple, Union, Dict

import pandas as pd

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
    URL_NUCCORE = 'https://www.ncbi.nlm.nih.gov/nuccore/{id}'
    URL_PUBMED = 'https://pubmed.ncbi.nlm.nih.gov/{id}'
    MATCH_COLORS = {0: None, 1: 'grey', 2: 'lightgreen', 3: 'green'}

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the tool.
        :param camel: CAMEL instance
        """
        super().__init__('ResFinder Reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid.
        :return: None
        """
        if 'TSV_genes' not in self._tool_inputs:
            raise InvalidInputSpecificationError('ResFinder input (TSV_genes) is required.')
        if 'TSV_pheno_general' not in self._tool_inputs:
            raise InvalidInputSpecificationError('ResFinder phenotype input (TSV_pheno_general) is required.')
        if 'resfinder' not in self._input_informs:
            raise InvalidInputSpecificationError('ResFinder informs are required.')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection(ResFinderReporter.TITLE, subtitle=self._input_informs['resfinder']['_name'])
        header, data = self.__parse_input_file()
        data = self.__add_links_and_format_data(data)
        self.__add_pheno_table(section)
        self.__add_output_table(section, header, data)
        self.__add_explanation_matches(section)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]

    def __parse_input_file(self) -> Tuple[Dict[str, List[str]], Dict[str, List[List[str]]]]:
        """
        Parses the input file.
        :return: Input file header, input file data
        """
        header, data = {}, {}
        for key in ['TSV_genes', 'TSV_pheno_general', 'TSV_point', 'TSV_pheno_species']:
            if key not in self._tool_inputs:
                continue
            input_file = self._tool_inputs[key][0]
            with open(input_file.path) as handle:
                header[key] = handle.readline().strip().split('\t')
                data[key] = [line.strip().split('\t') for line in handle.readlines()]
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

    def __add_pheno_table(self, section: HtmlReportSection) -> None:
        """
        Adds the table with the phenotype information.
        :param section: Report section
        :return: None
        """

        for input_key in [k for k in self._tool_inputs if 'pheno' in k]:

            data_pheno = pd.read_table(self._tool_inputs[input_key][0].path, comment='#', names=[
                'Antimicrobial', 'Class', 'WGS-predicted phenotype', 'Match', 'Genetic background'])
            data_pheno = data_pheno[data_pheno['WGS-predicted phenotype'].apply(lambda x: not pd.isna(x))]

            # Replace NA by dashes
            data_pheno.fillna('-', inplace=True)

            # Sort by first two columns
            data_pheno.sort_values(by=['Class', 'Antimicrobial'], inplace=True)

            # Get the colors for the resistance column, then drop the 'Match' column
            row_colors = data_pheno['Match'].apply(lambda x: ResFinderReporter.MATCH_COLORS[int(x)])
            data_pheno.pop('Match')

            # Convert to list format
            data_table = data_pheno.values.tolist()

            # Update the phenotype by a colored cell
            for row, color in zip(data_table, row_colors):
                index_wgs_pheno = data_pheno.columns.get_loc('WGS-predicted phenotype')
                row[index_wgs_pheno] = HtmlTableCell(row[index_wgs_pheno], color=color)

            if 'general' in input_key:
                section.add_header('Predicted phenotype (general)', 3)
                div = HtmlExpandableDiv('phenotype_general', 'general')
                div.add_table(data_table, data_pheno.columns, [('class', 'data')])
                section.add_html_object(div)
            else:
                section.add_header('Predicted phenotype (species-specific)', 3)
                section.add_table(data_table, data_pheno.columns, [('class', 'data')])

        section.add_warning_message(
            "The phenotype 'No resistance' should be interpreted with caution, as it only means that nothing in the "
            "used database indicate resistance, but resistance could exist from 'unknown' or not yet implemented "
            "sources.")

        # Add download links
        for input_key, title in zip(['TSV_pheno_species', 'TSV_pheno_general'], ['species-specific', 'general']):
            if input_key not in self._tool_inputs:
                continue
            relative_path = Path('resfinder', self._tool_inputs[input_key][0].path.name)
            section.add_file(self._tool_inputs[input_key][0].path, relative_path)
            section.add_link_to_file(f'Download {title} overview (TSV)', relative_path)

    def __add_output_table(
            self, section: HtmlReportSection, header: Dict[str, List[str]],
            data: Dict[str, List[Union[str, HtmlTableCell]]]) -> None:
        """
        Adds the output table.
        :param section: Report section
        :param header: Output table header
        :param data: Output table data
        :return: None
        """
        key_to_info = {'TSV_genes': ['resfinder', 'Detected AMR genes'],
                       'TSV_point': ['pointfinder', 'Detected AMR-conferring mutations']}
        for key in data:
            table = data[key]
            if len(table) > 0:
                relative_path = Path('resfinder', self.__generate_output_filename(key_to_info[key][0]))
                section.add_header(f'{key_to_info[key][1]}', 4)
                section.add_table(table, header[key], [('class', 'data')])
                section.add_file(self._tool_inputs[key][0].path, relative_path)
                section.add_link_to_file('Download (TSV)', relative_path)
            else:
                section.add_header(f'{key_to_info[key][1]}', 4)
                section.add_paragraph('No genes found.')

    def __add_links_and_format_data(self, data: Dict[str, List[List[str]]]) -> Dict[str, List[Union[str, HtmlTableCell]]]:
        """
        Adds PubMed links for mutations that have an associated PMID.
        Also formats data for taking into account color and floating points.
        :param data: Data
        :return: Data with links added
        """
        table_results = {}
        for key in ['TSV_genes', 'TSV_point']:
            table = data[key]
            edited_data = []
            if table:
                n_fields = len(table[0])
                for i in range(len(table)):
                    row = table[i]
                    if row[n_fields - 1] is not None:
                        access = row[n_fields - 1]
                        if key == 'TSV_genes':
                            publi = HtmlTableCell(str(access), link=ResFinderReporter.URL_NUCCORE.format(id=access))
                        elif key == 'TSV_point':
                            publi = HtmlTableCell(str(access), link=ResFinderReporter.URL_PUBMED.format(id=access))
                        row[n_fields - 1] = publi
                    try:
                        row[3] = f'{eval(row[3]):.2f}'
                    except NameError:
                        pass
                    if key == 'TSV_genes':
                        if row[3] == '100.00':
                            if row[1] == '100.00':
                                row = [HtmlTableCell(k, color='green') for k in row]
                            else:
                                row = [HtmlTableCell(k, color='lightgreen') for k in row]
                        else:
                            row = [HtmlTableCell(k, color='grey') for k in row]
                    else:
                        row = [HtmlTableCell(k, color='green') for k in row]
                    edited_data.append(row)
                table_results[key] = edited_data
        return table_results

    def __add_explanation_matches(self, section: HtmlReportSection) -> None:
        """
        Adds information about the different type of matches to the bottom of the report.
        :param section: Report section
        :return: None
        """
        section.add_header('Extra information', 3)
        section.add_paragraph('The following colors are used to denote the different type of hits:')
        section.add_table([
            [HtmlTableCell('', color='green'), 'Perfect match (100% over full length)'],
            [HtmlTableCell('', color='lightgreen'), 'Coverage 100%, identity <100%'],
            [HtmlTableCell('', color='grey'), 'Coverage <100%', 'identity <= 100%'],
            [HtmlTableCell('', color=None), 'No match found'],
        ], None, [('class', 'data')])
