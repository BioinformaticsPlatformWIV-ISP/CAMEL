import logging
from typing import List

import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlexpandablediv import HtmlExpandableDiv
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class GenomicContext(Tool):
    """
    Tool to link detected genes to their genomic context (e.g. plasmid / chromosome).
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Genomic context', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'dbs' not in self._input_informs:
            raise InvalidInputSpecificationError('Database informs input is required')
        for data_db in self._input_informs['dbs']:
            if f"TSV_{data_db['key']}" not in self._tool_inputs:
                raise InvalidInputSpecificationError(f"TSV file for '{data_db['key']}' is required")
        super()._check_input()

    def _get_plasmid_status(self, contig_name: str, plasmids: List[str]) -> List[HtmlTableCell]:
        """
        Returns the plasmid status for the given contig.
        :param contig_name: Name of the contig
        :param plasmids: List of plasmids
        :return: List of cells for the output table
        """
        contig_status = self._input_informs['mob_recon']['contig_report'].get(contig_name)
        cells = [HtmlTableCell('X', color='green') if contig_status is None else HtmlTableCell('-')]
        for plasmid in plasmids:
            cells.append(HtmlTableCell('X', color='green') if plasmid == contig_status else HtmlTableCell('-'))
        return cells

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse mob_recon input data
        plasmids = sorted(list(x for x in set(
            self._input_informs['mob_recon']['contig_report'].values()) if x is not None))

        # Create the output report
        section = HtmlReportSection('Genomic context', level=3)

        # Check if the detection method is BLAST
        if self._parameters['detection_method'].value != 'blast':
            logging.warning('Genomic context can only be predicted with blast as detection method')
            section.add_paragraph('Predicting genomic context is only performed when the detection method is blast.')
        else:
            if self._parameters['read_type'] == 'illumina':
                section.add_warning_message(
                    'Predicting genomic context based solely on short-read data is error-prone and should only be '
                    'considered as an indication.')

            # Create rows
            for db in self._input_informs['dbs']:
                try:
                    data_db = pd.read_table(self._tool_inputs[f"TSV_{db['key']}"][0].path)
                except IndexError:
                    logging.warning(f"No hits found for database {db['key']}")
                    continue
                section.add_header(db['title'], 3)
                if len(data_db) > 10:
                    div = HtmlExpandableDiv(f"genomic_context-{db['key']}", f'{len(data_db)} rows.')
                else:
                    div = HtmlElement('div')
                div.add_table([
                    [f"<i>{row[db['gene']]}</i>",
                     *self._get_plasmid_status(row[db['contig']], plasmids)] for row in data_db.to_dict('records')
                ], ['Key', 'Chromosome', *plasmids], [('class', 'data')])
                section.add_html_object(div)

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
