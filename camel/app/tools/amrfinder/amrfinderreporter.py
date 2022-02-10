from pathlib import Path
from typing import List

import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class AMRFinderReporter(Tool):
    """
    Reporting class for the AMRFinder tool.
    """

    URL_ACCESSION_BASE = f'https://www.ncbi.nlm.nih.gov/protein/'
    OUTPUT_COLS_OVERVIEW = [
        'Gene symbol', 'Element type', 'Element subtype', 'Class', 'Subclass', 'Method',
        'Accession of closest sequence']
    OUTPUT_COLS_ALN = [
        'Gene symbol', 'Contig id', 'Start', 'Stop', 'Strand', 'Target length', 'Alignment length',
        'Reference sequence length', '% Coverage of reference sequence', '% Identity to reference sequence']
    OUTPUT_COL_ACCESSION = 'Accession of closest sequence'

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the tool.
        :param camel: CAMEL instance.
        """
        super().__init__('AMRFinder reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("TSV input is required")
        if 'amrfinder' not in self._input_informs:
            raise InvalidInputSpecificationError("AMRFinder informs are required")
        super()._check_input()

    def __add_stats_table(self, section: HtmlReportSection, df_in: pd.DataFrame, columns: List[str]) -> None:
        """
        Formats the column data.
        :param section: Report section
        :param df_in: Input dataframe
        :param columns: Selected columns
        :return: None
        """
        # Reformat data
        rows_out = []
        for index, row in df_in.iterrows():
            row_current = []
            for col in columns:
                if col.startswith('%'):
                    # XX.XX% notation for percentages
                    row_current.append(HtmlTableCell(f'{row[col]:.2f}', color=row['row_color']))
                elif col == AMRFinderReporter.OUTPUT_COL_ACCESSION:
                    # Links for accession numbers
                    row_current.append(HtmlTableCell(
                        row[col], link=f'{AMRFinderReporter.URL_ACCESSION_BASE}{row[col]}', color=row['row_color']))
                else:
                    row_current.append(HtmlTableCell(row[col], color=row['row_color']))
            rows_out.append(row_current)

        # Add table
        section.add_table(rows_out, columns, [('class', 'data')])

    @staticmethod
    def __determine_row_color(row_data: pd.Series) -> str:
        """
        Determines the row color.
        :param row_data: Row data
        :return: Color
        """
        if all(row_data[x] == 100.0 for x in ('% Identity to reference sequence', '% Coverage of reference sequence')):
            return 'green'
        elif row_data['% Coverage of reference sequence'] == 100.0:
            return 'lightgreen'
        else:
            return 'grey'

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('AMRFinder', subtitle=self._input_informs['amrfinder']['_name'])

        # Parse input data
        data_out = pd.read_table(self._tool_inputs['TSV'][0].path)

        # Overview table
        if len(data_out) > 0:
            data_out['row_color'] = data_out.apply(lambda x: AMRFinderReporter.__determine_row_color(x), axis=1)
            section.add_header('Overview', level=3)
            data_out['Class'] = data_out['Class'].apply(lambda x: ', '.join([p.title() for p in x.split('/')]))
            data_out['Subclass'] = data_out['Subclass'].apply(lambda x: ', '.join([p.title() for p in x.split('/')]))
            self.__add_stats_table(section, data_out, AMRFinderReporter.OUTPUT_COLS_OVERVIEW)

            # Alignment stats table
            section.add_header('Alignment', level=3)
            self.__add_stats_table(section, data_out, AMRFinderReporter.OUTPUT_COLS_ALN)
        else:
            section.add_paragraph('No hits found')

        # Output and info section
        section.add_header('Output and extra information', level=3)
        relative_path = Path('amrfinder') / self._tool_inputs['TSV'][0].path.name
        section.add_link_to_file('Download complete results (TSV)', relative_path)
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)
        section.add_paragraph(f"Database version: {self._input_informs['amrfinder']['db_version']}")

        # Set the tool outputs
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
