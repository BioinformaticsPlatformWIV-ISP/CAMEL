
from pathlib import Path

import numpy as np
import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class GenotyphiReporter(Tool):
    """
    Parses Genotyphi csv output reports and generates a HTML report section for the final report.
    """

    TITLE = 'Genotyphi'

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Genotyphi Reporter', '0.1', camel)
        self._section = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._section = HtmlReportSection(GenotyphiReporter.TITLE,
                                          subtitle=self._input_informs['genotyphi']['_name'])
        self.__add_antibiotic_sensitivity()
        self._section.add_horizontal_line()
        self.__add_lineage_information()
        self._section.add_horizontal_line()
        self.__add_output_table_link()
        relative_path = Path('genotyphi', 'summary_out.tsv')
        self._section.add_file(self._tool_inputs['TSV_output'][0].path, relative_path)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]
        self._section.add_html_object(
            HtmlElement('a', 'Description of the fields of the tables',
                        [('href', 'https://github.com/Mykrobe-tools/mykrobe/wiki/AMR-prediction-output')]))
        self.__add_database_information()

    def __add_antibiotic_sensitivity(self) -> None:
        """
        Adds the table with the antibiotic sensitivity.
        :return: None
        """
        header = ['Antibiotic', 'Susceptibility', 'Variant', 'Genes']
        data = []
        results = pd.read_csv(self._tool_inputs['CSV'][0].path)
        # replace all nan by dashes
        results.replace(np.nan, '-', inplace=True)
        for i in range(results.shape[0]):
            data.append([results.iloc[i, 1], results.iloc[i, 2], results.iloc[i, 3], results.iloc[i, 4]])

        # start writing in the report the table and the headers
        self._section.add_header('Antibiotic susceptibility', 3)
        self._section.add_table(data, header, [('class', 'data')])

    def __add_lineage_information(self) -> None:
        """
        Adds the table with the lineage information.
        :return: None
        """
        header = ["Phylo group", "Species", "Lineage", "Phylo group per covg", "Species per covg", "Lineage per covg",
                  "Phylo group depth", "Species depth", "Lineage depth"]
        results = pd.read_csv(self._tool_inputs['CSV'][0].path)
        # replace all nan by dashes
        results.replace(np.nan, '-', inplace=True)
        data = [results.iloc[0, 10:19]]
        # display the numbers with decimal as two decimal numbers
        if isinstance(data[0][2], float):
            data[0][2] = f"{data[0][2]:.2f}"
        else:
            data[0][2] = '-'
        data[0][3] = f"{data[0][3]:.2f}"
        # start writing in the report the table and the headers
        self._section.add_header('Lineage information for <i>Salmonella typhi</i>', 3)
        self._section.add_table(data, header, [('class', 'data')])

    def __add_output_table_link(self) -> None:
        """
        Adds link to the output table (tsv) for this assay.
        :return: None
        """
        relative_path = Path('genotyphi', 'summary_out.tsv')
        self._section.add_link_to_file("Download (TSV)", relative_path)

    def __add_database_information(self) -> None:
        """
        Adds the database information to the report.
        :return: None
        """
        self._section.add_paragraph('Last updated: {}'.format(self._input_informs['genotyphi'].get(
            'last_update_date', '{LAST_UPDATE_DATE}')))