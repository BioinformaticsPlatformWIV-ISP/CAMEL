import pandas as pd

from camel.app.camel import Camel
from camel.app.tools.tool import Tool
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.io.tooliovalue import ToolIOValue


class CharacterizeNeisseriaCapsuleReporter(Tool):
    """
    Creates an HTML output report for the characterize_neisseria_capsule tool.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the reporter.
        """
        super().__init__('CharacterizeNeisseriaCapsule reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Serogrouping tool tabular output is required ('TSV')")
        super()._check_input()

    def __add_serogroup_prediction(self, section: HtmlReportSection) -> None:
        """
        Adds the detected serogroup to the output report section.
        :param section: Report output section
        :return: None
        """
        # Parse the input file
        table_data = pd.read_table(self._tool_inputs['TSV'][0].path)

        # Extract serogroup
        assigned_SG = table_data['SG'].iloc[0]

        # Create table data
        section.add_table(
            [[HtmlTableCell(assigned_SG, CharacterizeNeisseriaCapsuleReporter.__get_color(assigned_SG))]],
            ['Predicted serogroup'], [('class', 'data')])

    def __add_serogroup_table(self, section: HtmlReportSection) -> None:
        """
        Adds a table with the detected serogroup to the output report section.
        :param section: Report output section
        :return: None
        """
        # Parse the input file
        data = pd.read_table(self._tool_inputs['TSV'][0].path)

        # Remove the first column (contains sample ID)
        data.pop('Query')

        # Extract note to display it under the table
        note_txt = data['Notes'].iloc[0]
        data.pop('Notes')

        # Rename columns
        header = ['Predicted serogroup', 'Identified capsule genes']

        # Create table data
        table_data = []
        for values in data.itertuples(index=False, name=None):
            row = list(values)

            # Color the cell with the predicted serogroup
            row[0] = HtmlTableCell(row[0], CharacterizeNeisseriaCapsuleReporter.__get_color(row[0]))
            table_data.append(row)
        section.add_table(table_data, header, [('class', 'data')])
        section.add_paragraph(f'<b>Note:</b> {note_txt}.')


    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Create overview table with status
        section = HtmlReportSection('Characterize Neisseria Capsule tool')
        self.__add_serogroup_prediction(section)

        # Create table with hits
        section.add_header('Serogroup information', 3)
        self.__add_serogroup_table(section)

        # Store the output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    @staticmethod
    def __get_color(serogroup: str) -> str:
        """
        Returns the color for the assigned serogroup.
        :param serogroup: tool predicted serogroup
        :return: Serogroup cell
        """
        if serogroup == 'NG':
            return 'yellow'
        return 'green'
