import re
from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class MOBReconReporter(Tool):
    """
    Creates an HTML report for MOB-recon.
    """

    COLUMN_MAPPING = {
        'id': {'title': 'ID'},
        'num_contigs': {'title': 'Nb. of contigs'},
        'size': {'title': 'Size', 'fmt': lambda x: f'{x:,}'},
        'gc': {'title': '% GC-content', 'fmt': lambda x: f'{x * 100:.2f}'},
        'rep_type(s)': {'title': 'Rep. types', 'fmt': lambda x: x.replace(',', ', ')},
        'relaxase_type(s)': {'title': 'Relaxase types', 'fmt': lambda x: x.replace(',', ', ')}
    }

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('MOB-recon reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError('TSV input is required')
        if 'TSV_contigs' not in self._tool_inputs:
            raise InvalidInputSpecificationError('TSV contigs input is required')
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTA input is required')
        if 'mob_recon' not in self._input_informs:
            raise InvalidInputSpecificationError("MOB-recon informs input is required")
        super()._check_input()

    def _add_overview_table(self, section: HtmlReportSection) -> None:
        """
        Adds an overview table to the report.
        :param section: Report section
        :return: None
        """
        if len(self._input_informs['mob_recon']['detected_plasmids']) > 0:
            data_overview = pd.read_table(self._tool_inputs['TSV'][0].path)
        else:
            section.add_paragraph('No plasmids detected')
            return

        data_overview['id'] = data_overview['sample_id'].apply(lambda x: re.search('.*:(.*)', x).group(1))
        section.add_header('Overview', 3)
        table_data = [
            [d.get('fmt', lambda x: x)(row[col]) for
             col, d in MOBReconReporter.COLUMN_MAPPING.items()] for row in data_overview.to_dict('records')
        ]

        # Add column with download link for FASTA files
        for row in table_data:
            id_ = row[0]
            path_fasta = next(io.path for io in self._tool_inputs['FASTA'] if id_ in io.path.name)
            relative_path = Path('mob-suite', path_fasta.name)
            section.add_file(path_fasta, relative_path)
            row.append(HtmlTableCell('Download (FASTA)', link=str(relative_path)))

        # Add table
        section.add_table(
            table_data, [c['title'] for c in MOBReconReporter.COLUMN_MAPPING.values()] + ['Sequence'],
            [('class', 'data')])

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('MOB-recon', subtitle=self._input_informs['mob_recon']['_name'])
        self._add_overview_table(section)

        # Download overview
        relative_path = Path('mob-suite', 'mob_recon.tsv')
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)
        section.add_link_to_file('Download complete output (TSV)', relative_path)

        # Download contig report
        relative_path = Path('mob-suite', 'contig_report.tsv')
        section.add_file(self._tool_inputs['TSV_contigs'][0].path, relative_path)
        section.add_link_to_file('Download contig report (TSV)', relative_path)

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
