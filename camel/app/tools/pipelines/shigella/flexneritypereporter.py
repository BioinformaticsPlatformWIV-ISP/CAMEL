from pathlib import Path
from typing import List, Any

from camel.app.camel import Camel
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlexpandablediv import HtmlExpandableDiv
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.pipelines.shigella.flexneritypedetector import COORD_MINUS35_BOX, COORD_MINUS10_TA_BOX
from camel.app.tools.tool import Tool


class FlexneriTypeReporter(Tool):
    """
    This tool is used to create a report for the Shigella Flexneri type determination.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance.
        :return: None
        """
        super().__init__('Shigella: flexneri type reporter', '0.1', camel)
        self._section = HtmlReportSection('Flexneri type determination')
        self._sub_folder = Path('subspecies_identification')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._section.add_paragraph('Detected flexneri type: <b>{}</b>'.format(
            self._input_informs['detection']['detected_type']))
        if self._input_informs['subspecies']['detected_subspecies'] != 'flexneri':
            self._section.add_warning_message('Subspecies might not be <i>flexneri</i>.')
        self.__add_profiles_table()

        self._section.add_header('Hits', 4)
        table_data = []
        for locus, info in sorted(self._input_informs['detection']['loci'].items()):
            table_data.append([
                locus,
                HtmlTableCell('Yes', color='green') if info['detected'] is True else HtmlTableCell('No', color='red'),
                ', '.join(info['mutations']['stop']) if len(info['mutations']['stop']) > 0 else '-',
                ', '.join(info['mutations']['frameshift']) if len(info['mutations']['frameshift']) > 0 else '-',
                HtmlTableCell('Download (VCF)', link=str(self._section.add_file(
                    info['VCF'], self._sub_folder / f'variants-csq_{locus}.vcf'))) if info['detected'] is True else '-'
            ])
        header = ['Locus', 'Detected', 'Stop mutations', 'Frameshift mutations', 'VCF']
        self._section.add_table(table_data, header, [('class', 'data')])

        self._section.add_header('<i>gtr</i> promotor', 4)
        table_data = [
            ['Wild type <i>gtr</i> promotor:', self._input_informs['detection']['wt_gtr_promotor']],
            ['<i>gtr</i> promotor depth:', '{}X'.format(self._input_informs['gtr_depth']['median_depth'])]
        ]
        self._section.add_table(table_data, table_attributes=[('class', 'information')])

        self._section.add_paragraph(
            '''The promotor is considered wild type (WT) if there are no mutations in the -35 box or in the -10 TA box.
            If there is any mutation in one of those regions the <i>gtrX</i> gene is considered not to be expressed 
            for the determination of the <i>flexneri</i> type.
            ''')

        # Add table with mutations
        header = ['Position', 'Ref', 'Alt']
        table_data = []
        for key, title in (
                ('-35_box', '-35 box ({} -> {})'.format(*COORD_MINUS35_BOX)),
                ('-10_TA_box', '-10 TA box ({} -> {})'.format(*COORD_MINUS10_TA_BOX)),
                ('other', 'Other')):
            table_data.extend(self.__generate_promotor_section_table_rows(key, title, len(header)))
        div_promotor_mutations = HtmlExpandableDiv('promotor_mutations', 'promotor mutations')
        div_promotor_mutations.add_table(table_data, header, [('class', 'data')])
        self._section.add_html_object(div_promotor_mutations)

        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __generate_promotor_section_table_rows(self, key: str, title: str, nb_of_columns: int) -> List[List[Any]]:
        """
        Generates the section for the table with the given mutation.
        :param key: Key in informs
        :param title: Title in table
        :param nb_of_columns: Number of columns
        :return: List of table rows
        """
        rows = [[HtmlElement('th', title, attributes=[('colspan', nb_of_columns)])]]
        if len(self._input_informs['detection']['promotor_variants'][key]) > 0:
            for relative_pos, variant in self._input_informs['detection']['promotor_variants'][key]:
                rows.append([relative_pos, variant.REF, str(variant.ALT[0])])
        else:
            rows.append([HtmlElement('td', 'No mutations found', attributes=[('colspan', nb_of_columns)])])
        return rows

    def __add_profiles_table(self) -> None:
        """
        Adds the table with the Flexneri type profiles.
        :return: None
        """
        self._section.add_header('Profiles', 4)
        div = HtmlExpandableDiv('flexneri_profiles', 'profiles')
        table_data = []
        for profile, loci in self._input_informs['detection']['profiles'].items():
            color = 'yellow' if self._input_informs['detection']['detected_type'] == profile else None
            row = [profile] + ['+' if present else '-' for _, present in loci.items()]
            table_data.append([HtmlTableCell(x, color) for x in row])
        header = ['Profile'] + list(self._input_informs['detection']['loci'].keys())
        div.add_table(table_data, header, [('class', 'data')])
        self._section.add_html_object(div)
