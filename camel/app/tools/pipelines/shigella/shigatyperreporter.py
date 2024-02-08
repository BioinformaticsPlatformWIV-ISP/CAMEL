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

        # Create table data
        hits_table = []
        for values in gene_hits.itertuples(index=False, name=None):
            row = list(values)
            hits_table.append(row)

        # Rename columns
        header = ['Hit', 'Number of reads', 'Length Covered', 'Reference length', '% covered',
                  'Number of variants', '% accuracy']

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
