import json
from typing import List, Optional

import os

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlexpandablediv import HtmlExpandableDiv
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class HtmlReporterTyping(Tool):
    """
    Tool that creates HTML reports for the sequence typing pipeline.

    Input:
        - HTML: Path to the HTML file to write the report
        - DIR: Directory to store files that are included in the HTML report
        - Informs 'Scheme': Information about the scheme
        - VAL_Hits: Hits detected for each locus
    Output:
        - HTML: Path to the generated report
    """

    INFO_FILENAME = 'sequence_typing.json'

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('HTML Reporter', '0.1', camel)
        self._report_section = None
        self._output_folder = None
        self._sub_folder = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__initialize_report()
        if 'ST' in self._input_informs:
            self.__add_sequence_type()

        # Add tables with the detected hits
        add_subtitle = True if all(len(self._tool_inputs[f'hits_{key}']) > 1 for key in ('nucl', 'pept')) else False
        if len(self._tool_inputs['hits_nucl']) != 0:
            self.__add_output_table(self._tool_inputs['TSV_nucl'][0].path, self._tool_inputs['hits_nucl'],
                                    'Nucleotide loci' if add_subtitle else None)
        if len(self._tool_inputs['hits_pept']) != 0:
            self.__add_output_table(self._tool_inputs['TSV_pept'][0].path, self._tool_inputs['hits_pept'],
                                    'Peptide loci' if add_subtitle else None)
        if 'forced_detection_method' in self._parameters:
            self._report_section.add_alert(
                f"Allele detection performed with <b>{self._parameters['forced_detection_method'].value}</b>.", 'info')
        self._report_section.add_paragraph('Last updated: {}'.format(self._input_informs['scheme']['last_updated']))
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section)]
        self.__export_analysis_metadata()

    def __initialize_report(self) -> None:
        """
        Initializes the HTML report.
        :return: None
        """
        self._report_section = HtmlReportSection(self._input_informs['scheme']['title'], 3)
        self._sub_folder = os.path.join(
            'sequence_typing', FileSystemHelper.make_valid(self._input_informs['scheme']['name']))

    def __add_sequence_type(self) -> None:
        """
        Adds the sequence type to the report.
        :return: None
        """
        profile = self._input_informs['ST']['sequence_type']
        header = [key for key, _ in profile.metadata]
        table_data = [[value if value != '' else '-' for _, value in profile.metadata]]
        st = table_data[0][0]
        table_data[0][0] = HtmlTableCell(st, 'green' if st != '-' else 'red')
        self._report_section.add_table(table_data, header, table_attributes=[('class', 'data')])

    def __add_output_table(self, output_tsv: str, hits_io: List[ToolIOValue], sub_header: Optional[str]) -> None:
        """
        Adds the output table with the detected alleles.
        :param output_tsv: Tabular output file
        :param hits_io: Detected hits
        :param sub_header: If not None, this sub header is added to the report
        :return: None
        """
        table_header = hits_io[0].value.html_column_names()
        table_data = [h.value.to_html_row(self._report_section, self._sub_folder) for h in sorted(
            hits_io, key=lambda x: x.value.locus)]

        if sub_header is not None:
            self._report_section.add_header(sub_header, 4)

        # Add slider for big tables
        if len(hits_io) > 12:
            div = HtmlExpandableDiv('table-{}'.format(
                self._input_informs['scheme']['name'].lower()), f'alleles ({len(hits_io)})')
            div.add_table(table_data, table_header, [('class', 'data')])
            self._report_section.add_html_object(div)
        else:
            self._report_section.add_table(table_data, table_header, [('class', 'data')])

        relative_path = os.path.join(self._sub_folder, os.path.basename(output_tsv))
        self._report_section.add_file(output_tsv, relative_path)
        self._report_section.add_link_to_file("Download (TSV)", relative_path)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'scheme' not in self._input_informs:
            raise InvalidInputSpecificationError("Scheme information is required")
        super()._check_input()

    def __export_analysis_metadata(self) -> None:
        """
        Exports the analysis metadata file. The information can be used for further processing of the sequence typing
        output (e.g. for generating MLST trees).
        :return: None
        """
        path = os.path.join(self._folder, HtmlReporterTyping.INFO_FILENAME)
        with open(path, 'w') as handle:
            json.dump({
                'scheme': self._input_informs['scheme']['name'],
                'sample': self._tool_inputs['VAL_SAMPLE'][0].value
            }, handle)
        self._report_section.add_file(path, os.path.join(self._sub_folder, HtmlReporterTyping.INFO_FILENAME))
