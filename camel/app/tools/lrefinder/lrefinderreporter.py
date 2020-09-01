from typing import Union

import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class LREFinderReporter(Tool):
    """
    Tool to generate HTML output for the LRE-Finder tool.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('LRE-Finder reporter', '1.0', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'lrefinder' not in self._input_informs:
            raise InvalidInputSpecificationError('LRE-Finder informs are required')
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
        section = HtmlReportSection('LRE-Finder', subtitle=self._input_informs['lrefinder']['_name'])
        informs = self._input_informs['lrefinder']
        section.add_paragraph(f"Detected species: <i>{informs['species']}</i>")
        data_genes = pd.DataFrame(informs['genes'])
        section.add_header('Detected genes', 3)
        section.add_table([
            [LREFinderReporter.__format_value(row[col]) for col in data_genes.columns] for
            row in data_genes.to_dict('records')], data_genes.columns, [('class', 'data')])
        section.add_header('Detected mutations', 3)
        data_mutations = pd.DataFrame(informs['mutations'])
        section.add_table([
            [LREFinderReporter.__format_value(row[col]) for col in data_mutations.columns] for
            row in data_mutations.to_dict('records')], data_mutations.columns, [('class', 'data')])
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
