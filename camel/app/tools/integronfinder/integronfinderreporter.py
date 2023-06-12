from pathlib import Path
from typing import Union

import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class IntegronFinderReporter(Tool):
    """
    Tool to generate HTML output for the IntegronFinder tool.
    """

    COLUMNS = [
        {'key': 'ID_integron', 'title': 'Integron'},
        {'key': 'ID_replicon', 'title': 'Replicon'},
        {'key': 'pos_beg', 'title': 'Start'},
        {'key': 'pos_end', 'title': 'End'},
        {'key': 'strand', 'title': 'Strand'},
        {'key': 'type_elt', 'title': 'Element'},
        {'key': 'annotation', 'title': 'Annotation'},
        {'key': 'model', 'title': 'Model'},
        {'key': 'type', 'title': 'Type'},
        {'key': 'considered_topology', 'title': 'Topology'}
    ]

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('IntegronFinder reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'integron_finder' not in self._input_informs:
            raise InvalidInputSpecificationError('IntegronFinder informs are required')
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError('IntegronFinder output file is required (TSV)')
        super()._check_input()

    @staticmethod
    def __format_value(value: Union[str, float]) -> str:
        """
        Formats a value for the output table.
        :param value: Input value
        :return: Formatted value
        """
        if type(value) is float:
            return f'{value:.2f}'
        return value

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('IntegronFinder', subtitle=self._input_informs['integron_finder']['_name'])

        # Parse input data
        data_integrons = pd.read_table(self._tool_inputs['TSV'][0].path, comment='#')
        data_integrons.fillna('-')
        data_integrons['strand'] = data_integrons['strand'].apply(lambda x: '-' if x == -1 else '+')

        # Add output table
        section.add_table([
            [row[col['key']] for col in IntegronFinderReporter.COLUMNS] for
            row in data_integrons.to_dict('records')],
            [c['title'] for c in IntegronFinderReporter.COLUMNS],
            [('class', 'data')])

        # Add download link
        relative_path = Path('integron_finder', f"integrons_{self._parameters['name'].value}.tsv")
        section.add_link_to_file('Download (TSV)', relative_path=relative_path)
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
