import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class SpeciesDeterminationReporter(Tool):
    """
    This class is used to generate an output report for the species determination.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Yersinia: species determination reporter', '1.0', camel)
        self._section = HtmlReportSection('Species determination')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        if self._input_informs['analysis']['found_match']:
            table_df = pd.read_csv(self._tool_inputs['TSV_analysis'][0].path, sep="\t", header=0).fillna('NA')
            species, lineage, biotype, serotype, match, threshold = self._input_informs['analysis']['best_match'].values()
            self._section.add_paragraph(f'<b>Species:</b> {species}')
            self._section.add_paragraph(f'<b>Lineage:</b> {lineage}')
            self._section.add_paragraph(f'<b>Biotype:</b> {biotype}')
            self._section.add_paragraph(f'<b>Serotype:</b> {serotype}')
            self._section.add_paragraph('All matches:')
            self.__add_table_detected_species(table_df)
        else:
            self._section.add_paragraph('No Yersinia species found.')
        self._section.add_paragraph("Species and lineage designations as defined by <a href=\"https://doi.org/10.1099%2Fmgen.0.000301\">Savin et al</a>.")
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_table_detected_species(self, table_df) -> None:
        """
        Adds the table with the detected species.
        :return: None
        """
        table_html = []
        #TODO: color code?
        for _, r in table_df.iterrows():
            match = r['match']
            threshold = r['threshold']
            row = [r['species'],r['lineage'],r['biotype'],r['serotype'], f'{match:.2f}', f'{threshold:.2f}']
            table_html.append(row)
        header = ['Species', 'Lineage', 'Biotype', 'Serotype', 'Proportion of loci matched', 'Threshold for species and lineage']
        self._section.add_table(table_html, header, [('class', 'data')])
