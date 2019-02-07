import json

import os

from camel.app.command.command import Command
from camel.app.components.html.htmlexpandablediv import HtmlExpandableDiv
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class SpoTypingReporter(Tool):
    """
    This tool creates a report for the SpoTyping output.
    """

    TITLE = 'Spoligotyping'

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL
        """
        super().__init__('SpoTyping Reporter', '0.1', camel)
        self._section = None
        self._section_spacers = HtmlReportSection('Spacer counts')

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'VAL_type_binary' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Binary spoligotype input (VAL_type_binary) is required")
        if 'VAL_type_octal' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Octal spoligotype input (VAL_type_octal) is required")
        if 'metadata' not in self._input_informs:
            raise InvalidInputSpecificationError("Spoligotype metadata is required")
        super()._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self._section = HtmlReportSection(SpoTypingReporter.TITLE, subtitle=self._input_informs['metadata']['version'])
        metadata = self.__extract_metadata(self._tool_inputs['VAL_type_octal'][0].value)
        table_data = [
            ['Binary:', self._tool_inputs['VAL_type_binary'][0].value],
            ['Octal:', self._tool_inputs['VAL_type_octal'][0].value],
            ['SIT number:', metadata['SIT']],
            ['SpolDB4 occurence:', metadata['total']],
            ['Label:', metadata['label']]
        ]
        self._section.add_table(table_data, table_attributes=[('class', 'information')])
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

        if 'LOG' in self._tool_inputs:
            self.__add_spacer_section()
        self.__add_parameter_section()

    def __add_spacer_section(self):
        """
        Creates a section that contains the number of times each spacer appeared.
        :return: None
        """
        header = ['Spacer', 'Error free hits', '1-error tolerant hits', 'Detected']
        with open(self._tool_inputs['LOG'][0].path) as handle:
            table_data = []
            for line in handle.readlines()[3:]:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                code = parts[3]
                detected = 'Yes' if code == '1' else 'No'
                color = 'green' if code == '1' else 'red'
                table_data.append([parts[0], parts[1], parts[2], HtmlTableCell(detected, color)])
        div = HtmlExpandableDiv('table-spacers', 'spacer counts')
        div.add_table(table_data, header, [('class', 'data')])
        self._section.add_html_object(div)

    def __extract_metadata(self, type_octal):
        """
        Extracts the metadata for the detected Spoligotype.
        :return: Spoligotype metadata
        """
        command = Command('{} echo $SPOTYPING_METADATA'.format(self._build_dependencies()))
        command.run_command(self._folder)
        metadata_path = command.stdout.strip()
        if not os.path.isfile(metadata_path):
            raise FileNotFoundError("No spoligotype metadata found")
        with open(metadata_path) as handle:
            metadata = json.load(handle)
        keys = ('SIT', 'geo', 'label', 'total')
        if type_octal in metadata:
            return {k: metadata[type_octal][k] for k in keys}
        else:
            return {k: 'NA' for k in keys}

    @staticmethod
    def generate_empty_section() -> HtmlReportSection:
        """
        Returns a report that is used when this analysis is disabled.
        :return: Report section
        """
        section = HtmlReportSection(SpoTypingReporter.TITLE)
        section.add_paragraph('Analysis disabled')
        return section

    def __add_parameter_section(self) -> None:
        """
        Adds a section with the parameters used for SpoTyping.
        :return: None
        """
        self._section.add_header('Parameters', 4)
        table_data = [
            ['Min. number strict hits', self._input_informs['spoligo_param']['min_strict']],
            ['Min. number relaxed hits', self._input_informs['spoligo_param']['min_relaxed']],
            ['Downsample factor', self._input_informs['spoligo_param']['downsample_factor']]
        ]
        self._section.add_table(table_data, ['Name', 'Value'], [('class', 'data')])
