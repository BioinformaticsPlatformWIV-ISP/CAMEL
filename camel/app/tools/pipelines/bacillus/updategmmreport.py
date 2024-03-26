import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class UpdateGMMReport(Tool):
    """
    Updates the GMM gene detection report with warning about detected GMMs.
    """
    INPUT_KEYS = ['TSV_STRAINS', 'TSV_GMM', 'VAL_HTML', 'TSV_GMM_DB']
    color_code = {'STRAIN_MATCH': 'green', 'GMM_MATCH': 'yellow', 'BOTH_MATCH': 'red'}

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the tool.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('UpdateGMMReport', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if any(key not in self._tool_inputs for key in self.INPUT_KEYS):
            raise InvalidInputSpecificationError(
                "Tool requires {} inputs".format(', '.join(UpdateGMMReport.INPUT_KEYS)))
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the tool.
        :return: None
        """
        self._parse_tsv_files()
        output_report = self._update_report()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(output_report)]

    def _update_report(self) -> HtmlReportSection:
        """
        Updates the report with the GMM warning table.
        :return: Updated report section
        """
        matches = self._parse_tsv_files()
        current_report_section = self._tool_inputs['VAL_HTML'][0].value
        current_report_section.add_header('Interpretation', level=4)
        if not (matches['strain'] and matches['construct']):
            current_report_section.add_paragraph('No GMM construct detected')
            return current_report_section

        table_to_add = list(zip(matches['strain'], matches['construct']))
        column_names = ['strain', 'construct']
        current_report_section.add_table(table_to_add, column_names, [('class', 'data')])

        if matches['strain']:
            current_report_section.add_paragraph(f'The strain matches closely to <b>{matches["strain"][0]}</b> '
                                                 f'which has been used in GMMs.')
        else:
            current_report_section.add_paragraph('The strain does not match any known GMM strains in the database.')
        if matches['construct']:
            current_report_section.add_paragraph(f'The <b>{matches["construct"][0]}</b> transgenic construct '
                                                 f'was detected in the strain.')
        else:
            current_report_section.add_paragraph('No transgenic constructs from the database were detected in the strain.')

        current_report_section.add_warning_message('The pipeline uses a targeted approach, which means that constructs '
                                                   'and/or strains that are not in the database will be missed.')

        return current_report_section

    def _parse_tsv_files(self) -> dict:
        """
        Parses the TSV files passed as input.
        :return: Dictionary with match, or False if no match is found
        """
        tsv_gmm_db = pd.read_csv(self._tool_inputs['TSV_GMM_DB'][0].path)

        strain_hits = []
        for f in self._tool_inputs['TSV_STRAINS']:
            with open(f.path) as handle:
                for line in handle:
                    spl = line.strip().split()
                    if 'closest_strain' in spl[0]:
                        if spl[1] in tsv_gmm_db['strain'].tolist():
                            strain_hits.append(spl[1])

        gmm_hits = []
        with open(self._tool_inputs['TSV_GMM'][0].path) as handle:
            all_lines = handle.readlines()
            gmm_hits_list = eval(all_lines[0].strip().split('\t')[1])
            for entry in gmm_hits_list:
                if entry[1] in tsv_gmm_db['construct'].tolist():
                    gmm_hits.append(entry[1])

        return {'strain': strain_hits,
                'construct': gmm_hits}
