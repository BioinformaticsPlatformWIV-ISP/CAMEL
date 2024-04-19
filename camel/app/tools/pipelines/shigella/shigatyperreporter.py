from typing import Dict, Any

import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class ShigaTyperReporter(Tool):
    """
    Creates an HTML output report for the ShigaTyper tool.
    """

    COLS = [
        {'key': 'Hit'},
        {'key': 'Number of reads', 'fmt': lambda x: f'{x:,}'},
        {'key': 'Length Covered', 'fmt': lambda x: f'{x:,}'},
        {'key': 'reference length', 'fmt': lambda x: f'{x:,}'},
        {'key': '% covered', 'fmt': lambda x: f'{x:.2f}'},
        {'key': 'Number of variants', 'fmt': lambda x: f'{int(x):,}'},
        {'key': '% accuracy', 'fmt': lambda x: f'{x:.2f}'}
    ]

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the reporter.
        """
        super().__init__('ShigaTyper reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("ShigaTyper tabular output is required ('TSV')")
        if 'shigatyper' not in self._input_informs:
            raise InvalidInputSpecificationError("ShigaTyper informs are required")
        super()._check_input()

    def __add_main_output(self, section: HtmlReportSection) -> None:
        """
        Adds a table with the ShigaTyper output to the report section.
        :param section: Report output section
        :return: None
        """
        # Parse the input files
        main_output = pd.read_table(self._tool_inputs['TSV'][0].path)

        # Remove the sample ID column
        main_output.pop('sample')

        # Create table data
        main_table = []
        for values in main_output.itertuples(index=False, name=None):
            row = list(values)
            main_table.append(row)

        # Rename columns
        header = ['Prediction',	'ipaB', 'Notes']

        section.add_table(main_table, header, [('class', 'data')])

    @staticmethod
    def __format_cell(value: Any, col: Dict) -> HtmlTableCell:
        """
        Formats the corresponding table cell.
        :param value: Input value
        :param col: Column metadata
        :return: HTML table cell
        """
        if 'fmt' not in col or value == '-':
            return HtmlTableCell(str(value))
        return HtmlTableCell(col['fmt'](value))

    def __add_shigatyper_hits(self, section: HtmlReportSection) -> None:
        """
        Adds a table with the ShigaTyper hits details to the report section.
        :param section: Report output section
        :return: None
        """
        # Parse the input files
        gene_hits = pd.read_table(self._tool_inputs['TSV_HITS'][0].path)

        # Remove the sample ID column
        gene_hits.pop('Unnamed: 0')

        # Add potentially missing columns
        col_list = ['Hit', 'Number of reads', 'Length Covered', 'reference length',
                    '% covered', 'Number of variants', '% accuracy']
        missing_columns = [col for col in col_list if col not in gene_hits.columns]
        gene_hits[missing_columns] = '-'

        # Create table data
        hits_table = []
        for row in gene_hits.to_dict('records'):
            hits_table.append([ShigaTyperReporter.__format_cell(
                row[col['key']], col) for col in ShigaTyperReporter.COLS])

        # Rename columns
        header = [c['key'].title() for c in ShigaTyperReporter.COLS]
        section.add_table(hits_table, header, [('class', 'data')])

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Create overview table with status
        section = HtmlReportSection('ShigaTyper', subtitle=self._input_informs['shigatyper']['_name'])
        species = self._input_informs['shigatyper']['species']
        section.add_table(
            [[HtmlTableCell(species)]], ['Serotyping'], [('class', 'data')])

        # Create table with hits
        section.add_header('Overview', 3)
        self.__add_main_output(section)
        self.__add_shigatyper_hits(section)

        # Store the output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
