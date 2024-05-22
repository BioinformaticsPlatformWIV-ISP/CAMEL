import re
from pathlib import Path
from typing import Optional

import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class AbriTAMRReporter(Tool):
    """
    Parses abriTAMR report output and generates an HTML output report.
    """

    TITLE = 'abriTAMR'

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the tool.
        :param camel: camel instance
        """
        super().__init__('AbriTAMR Reporter', '0.1', camel)
        self._section = None

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        super(AbriTAMRReporter, self)._check_input()
        if 'VAL_SPECIES' not in self._tool_inputs:
            raise InvalidInputSpecificationError("VAL_SPECIES is required")
        elif not all(key in self._tool_inputs for key in ('TXT_MATCHES', 'TXT_PARTIALS')):
            raise InvalidInputSpecificationError("AbriTAMR run output files must be provided")

    def _execute_tool(self) -> None:
        """
        Executes the reporter.
        :rtype: None
        """
        self._section = HtmlReportSection(AbriTAMRReporter.TITLE,
                                          subtitle=self._input_informs['ABRITAMR_RUN']['_name'])
        self.__add_summaries_tables()
        if str(self._tool_inputs['VAL_SPECIES'][0]) == 'Salmonella':
            self.__add_antibiogram()
        else:
            self._section.add_header('Antibiogram', 3)
            self._section.add_paragraph(f"Not available for species '{self._tool_inputs['VAL_SPECIES'][0]}'")
        self.__add_output_table_link()
        self.__add_database_information()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_summaries_tables(self) -> None:
        """
        Adds summary tables for the antibiotic hit files not in the antibiogram.
        :return: None
        """
        self._section.add_header('Matches', 3)
        self._section.add_paragraph('>90% coverage & >90% identity')
        self.___add_summary_file_table(self._tool_inputs['TXT_MATCHES'][0].path, 'matches')
        self._section.add_paragraph('Genes annotated with * indicate >90% coverage and identity between 90% and 100%. '
                                    'No further annotation indicates that the gene recovered exhibits '
                                    '100% coverage and 100% identity to a gene in the gene catalog.')
        self._section.add_header('Partial matches', 3)
        self._section.add_paragraph('50-90% coverage & >90% identity')
        self.___add_summary_file_table(self._tool_inputs['TXT_PARTIALS'][0].path, 'partials')
        self._section.add_horizontal_line()

    def ___add_summary_file_table(self, summary_file_path: Path, matching_type: str) -> None:
        """
        Parses an AbriTAMR output summary file (TXT_PARTIALS or TXT_MATCHES) and
        adds the table to the section.
        :param summary_file_path: path of TXT_PARTIALS or TXT_MATCHES
        :param matching_type: str, matches or partials
        :return: None
        """
        header = ['Functional drug class', 'Genes']
        data = []

        # Read the summary file into a pandas DataFrame
        df = pd.read_table(summary_file_path)

        # Extract relevant data from DataFrame
        if len(df.columns) > 1:
            line1 = df.columns[1:]
            line2 = df.iloc[0, 1:].tolist()
            for functional_drug_class, gene in zip(line1, line2):
                data.append([functional_drug_class, gene])

        self._section.add_table(sorted(data), header, [('class', 'data')])
        relative_path_file = Path('abritamr', f"summary_{matching_type}.txt")
        self._section.add_file(summary_file_path, relative_path_file)
        self._section.add_link_to_file(f"Download {matching_type} (TSV)", relative_path_file)

    def __add_antibiogram(self) -> None:
        """
        Add the antibiogram to the report for Salmonella.
        :return: None
        """
        header = ['Antibiotic', 'Resistance mechanisms detected', 'Interpretation']
        data = []
        # looks like this:
        # abritamr_Ampicillin_ResMech     None detected
        # abritamr_Ampicillin_Interpretation      Susceptible
        # abritamr_Cefotaxime (ESBL)_ResMech      None detected
        # abritamr_Cefotaxime (ESBL)_Interpretation       Susceptible
        # ...
        results = pd.read_csv(self._tool_inputs['TSV_output'][0].path, delimiter='\t', header=None)
        # replace all nan by dashes
        results.fillna('-', inplace=True)
        for i in range(0, results.shape[0]-1, 2):
            interpretation = results.iloc[i+1, 1]
            color = self.___get_interpretation_color(interpretation)
            row = [re.sub("_ResMech|abritamr_", "", results.iloc[i, 0]),
                   results.iloc[i, 1],
                   results.iloc[i+1, 1]]
            row = [HtmlTableCell(x, color) for x in row]
            data.append(row)
        # start writing in the report the table and the headers
        self._section.add_header('Antibiogram for <i>Salmonella</i>', 3)
        self._section.add_table(data, header, [('class', 'data')])

    @staticmethod
    def ___get_interpretation_color(interpretation: str) -> Optional[str]:
        """
        For a given interpretation returns either green, red or None
        :param interpretation: Resistant or Susceptible or any other string
        :return: str ('green' or 'red') or None
        """
        if interpretation == 'Susceptible':
            color = 'green'
        elif interpretation == 'Resistant':
            color = 'red'
        else:
            color = None
        return color

    def __add_output_table_link(self) -> None:
        """
        Add link to the output table (tsv) for this assay.
        :return: None
        """
        relative_path = Path('abritamr', 'summary_out.tsv')
        self._section.add_file(self._tool_inputs['TSV_output'][0].path, relative_path)
        if str(self._tool_inputs['VAL_SPECIES'][0]) == 'Salmonella':
            relative_path = Path('abritamr', 'summary_out.tsv')
            self._section.add_link_to_file("Download (TSV)", relative_path)

    def __add_database_information(self) -> None:
        """
        Adds the database information to the report.
        :return: None
        """
        self._section.add_horizontal_line()
        self._section.add_header('Database info', level=4)
        self._section.add_paragraph('Last updated: {}'.format(self._input_informs['ABRITAMR_RUN'].get(
            'last_update_date', 'n/a')))
