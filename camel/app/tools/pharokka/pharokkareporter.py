from pathlib import Path
from typing import Dict, Any

import numpy as np
import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class PharokkaReporter(Tool):
    """
    Parses Pharokka output files and generates an HTML report section for the final report.
    """

    TITLE = 'Pharokka'

    COLS = [
        {'key': 'Contig'},
        {'key': 'Gene'},
        {'key': 'Hit'},
        {'key': 'Alignment score'},
        {'key': 'Sequence identity'},
        {'key': 'Start'},
        {'key': 'Stop'},
        {'key': 'Frame'}
    ]

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Pharokka Reporter', '0.1', camel)
        self._section = None

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'STATS' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Pharokka output is required ('pharokka_cds_functions.tsv')")
        if 'CARD' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Pharokka output is required ('CARD hits')")
        if 'VFDB' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Pharokka output is required ('VFDB hits')")
        if 'INPHARED' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Pharokka output is required ('pharokka_top_hits_mash_inphared.tsv')")
        if 'pharokka' not in self._input_informs:
            raise InvalidInputSpecificationError("Pharokka informs are required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection(PharokkaReporter.TITLE,
                                    subtitle=self._input_informs['pharokka']['_name'])

        # Summary metrics
        section.add_header('Annotation summary metrics', 3)
        self.__add_stats(section)

        # Additional information
        self.__add_output_table_link(section)

        # Antibiotic sensitivity
        if 'show_card' in self._parameters:
            section.add_header('Antibiotic susceptibility (CARD database)', 3)
            self.__add_antibiotic_sensitivity(section)

        # Virulence factors
        if 'show_vfdb' in self._parameters:
            section.add_header('Virulence factor identification (VFDB database)', 3)
            self.__add_virulence_factors(section)

        # Phage identification
        if 'inphared' in self._parameters:
            section.add_header('Phage identification (INPHARED database)', 3)
            self.__add_identification(section)

        # Store the output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    def __add_stats(self, section: HtmlReportSection) -> None:
        """
        Adds the table with the summary metrics.
        :return: None
        """
        data = pd.read_table(self._tool_inputs['STATS'][0].path)

        stats = []
        for values in data.itertuples(index=False, name=None):
            row = list(values)
            stats.append(row)

        # Rename columns
        header = ['Description', 'Count', 'Contig']
        section.add_table(stats, header, [('class', 'data')])

    def __add_antibiotic_sensitivity(self, section: HtmlReportSection) -> None:
        """
        Adds the table with CARD AMR hits.
        :return: None
        """
        data_hits = pd.read_table(self._tool_inputs['CARD'][0].path)

        if data_hits.empty:
            section.add_paragraph(self._input_informs['pharokka'].get('card_hits'))
        else:
            # Create table data
            hits_table = []
            for values in data_hits.itertuples(index=False, name=None):
                row = list(values)
                hits_table.append(row)

            # Rename columns
            header = [c['key'].title() for c in PharokkaReporter.COLS]
            section.add_table(hits_table, header, [('class', 'data')])

    def __add_virulence_factors(self, section: HtmlReportSection) -> None:
        """
        Adds the table with VFDB hits.
        :return: None
        """
        data_hits = pd.read_table(self._tool_inputs['VFDB'][0].path)

        if data_hits.empty:
            section.add_paragraph(self._input_informs['pharokka'].get('vfdb_hits'))
        else:
            # Create table data
            hits_table = []
            for values in data_hits.itertuples(index=False, name=None):
                row = list(values)
                hits_table.append(row)

            # Rename columns
            header = [c['key'].title() for c in PharokkaReporter.COLS]
            section.add_table(hits_table, header, [('class', 'data')])

    def __add_identification(self, section: HtmlReportSection) -> None:
        """
        Adds the table with the INPHARED information.
        :return: None
        """
        data = pd.read_table(self._tool_inputs['INPHARED'][0].path)

        # Create table data
        info_table = []
        for values in data.itertuples(index=False, name=None):
            row = list(values)
            info_table.append(row)

        # Rename columns
        header = ['Contig','Accession', 'mash_distance', 'mash_pval', 'mash_matching_hashes',
                  'Descriptio', 'Classification', 'Genome_Length_(bp)',
                  'Jumbophage', 'molGC_(%)', 'Molecule', 'Modification_Date', 'Number_CDS',
                  'Positive_Strand_(%)', 'Negative_Strand_(%)', 'Coding_Capacity_(%)',
                  'Low_Coding_Capacity_Warning', 'tRNAs', 'Host', 'Lowest_Taxa',
                  'Genus', 'Sub-family', 'Family', 'Order', 'Class', 'Phylum',
                  'Kingdom', 'Realm', 'Baltimore_Group', 'Genbank_Division',
                  'Isolation_Host_(beware_inconsistent_and_nonsense_values)']
        section.add_table(info_table, header, [('class', 'data')])

    def __add_output_table_link(self, section: HtmlReportSection) -> None:
        """
        Adds link to the genbank annotation file.
        :return: None
        """
        relative_path = Path('pharokka', 'pharokka.gbk')
        section.add_file(self._tool_inputs['GBK'][0].path, relative_path)
        section.add_link_to_file("Download (GBK)", relative_path)

    @staticmethod
    def __format_cell(value: Any, col: Dict[str, Any]) -> HtmlTableCell:
        """
        Formats the corresponding table cell.
        :param value: Input value
        :param col: Column metadata
        :return: HTML table cell
        """
        if 'fmt' not in col:
            return HtmlTableCell(str(value))
        return HtmlTableCell(col['fmt'](value))
