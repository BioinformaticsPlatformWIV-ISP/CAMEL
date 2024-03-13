from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html.htmlexpandablediv import HtmlExpandableDiv
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class ResFinderReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the ResFinder output.
    """

    TITLE = 'ResFinder'
    URL_NUCCORE = 'https://www.ncbi.nlm.nih.gov/nuccore/{id}'
    URL_PUBMED = 'https://pubmed.ncbi.nlm.nih.gov/{id}'
    MATCH_COLORS = {0: None, 1: 'grey', 2: 'lightgreen', 3: 'green'}

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the tool.
        :param camel: CAMEL instance
        """
        super().__init__('ResFinder Reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid.
        :return: None
        """
        if 'TSV_pheno_general' not in self._tool_inputs:
            raise InvalidInputSpecificationError('ResFinder phenotype input (TSV_pheno_general) is required.')
        if 'resfinder' not in self._input_informs:
            raise InvalidInputSpecificationError('ResFinder informs are required.')
        super()._check_input()

    def __add_phenotype_table(self, section: HtmlReportSection, key: str) -> None:
        """
        Adds the table with the phenotype overview.
        :param section: Report section
        :param key: Key (species, general)
        :return: None
        """
        # Parse input
        data_pheno = pd.read_table(
            self._tool_inputs[f'TSV_pheno_{key}'][0].path,
            names=['Antimicrobial', 'Class', 'WGS-predicted phenotype', 'Match', 'Genetic background'],
            comment='#'
        )
        data_pheno.sort_values(by=['Class', 'Antimicrobial'], inplace=True)

        # Add table to report
        overview_type = 'species-specific' if key == 'species' else 'general'
        header = ['Class', 'Antimicrobial', 'WGS-predicted phenotype', 'Genetic background']
        section.add_header(f'Predicted phenotype ({overview_type})', 3)
        div = HtmlExpandableDiv(f'pheno_{key}', 'overview')
        div.add_table([(
            row['Class'],
            row['Antimicrobial'],
            HtmlTableCell(row['WGS-predicted phenotype'], color=ResFinderReporter.MATCH_COLORS[row['Match']]),
            row['Genetic background'] if not pd.isna(row['Genetic background']) else '-'
        ) for row in data_pheno.to_dict('records')], header, [('class', 'data')])
        section.add_html_object(div)

        # Add download link
        relative_path = Path('resfinder', self._tool_inputs[f'TSV_pheno_{key}'][0].path.name)
        section.add_file(self._tool_inputs[f'TSV_pheno_{key}'][0].path, relative_path)
        section.add_link_to_file(f'Download {overview_type} overview (TSV)', relative_path)

    @staticmethod
    def __get_row_color(row: pd.Series) -> str:
        """
        Returns the color for the gene detection row.
        :param row: Input row
        :return: Color
        """
        if row['perc_cov'] == 100.0 and row['Identity'] == 100.0:
            return 'green'
        elif row['perc_cov'] == 100.0:
            return 'lightgreen'
        return 'grey'

    @staticmethod
    def __get_accession_cell(accession: str, color: str, is_pmid: bool = False) -> HtmlTableCell:
        """
        Returns a table cell with a link to the accession.
        :param accession: Accession numbers
        :param color: Cell color
        :param is_pmid: If true, the accession is a PubMed accession
        :return: Formatted table cell
        """
        if len(accession.strip()) == 0:
            return HtmlTableCell('-', color=color)
        url = (ResFinderReporter.URL_NUCCORE if not is_pmid else ResFinderReporter.URL_PUBMED).format(id=accession)
        return HtmlTableCell(str(accession), link=url, color=color)

    def __add_genes_table(self, section: HtmlReportSection) -> None:
        """
        Adds the table wit the detected AMR genes.
        :param section: Report section
        :return: None
        """
        # Parse input
        data_genes = pd.read_table(self._tool_inputs['TSV_genes'][0].path, na_values=['NA..NA'])
        logger.info(f'{len(data_genes)} genes parsed')
        cols_original = list(data_genes.columns)
        data_genes['perc_cov'] = data_genes['Alignment Length/Gene Length'].apply(
            lambda x: 100 * int(x.split('/')[0]) / int(x.split('/')[1]))
        data_genes['color'] = data_genes.apply(lambda row: ResFinderReporter.__get_row_color(row), axis=1)

        # Add table
        section.add_header('Detected AMR genes', 3)
        section.add_table([[
            *[HtmlTableCell(
                f'{row[col]:.2f}' if isinstance(row[col], float) else row[col],
                color=row['color']
            ) for col in cols_original if col != 'Accession no.'],
            ResFinderReporter.__get_accession_cell(row['Accession no.'], row['color'])
        ] for row in data_genes.fillna('-').to_dict('records')], cols_original, [('class', 'data')])

        # Add download link
        if len(data_genes) > 0:
            relative_path = Path('resfinder4', self._tool_inputs['TSV_genes'][0].path.name)
            section.add_file(self._tool_inputs['TSV_genes'][0].path, relative_path)
            section.add_link_to_file('Download (TSV)', relative_path)

    def __add_mutations_table(self, section: HtmlReportSection) -> None:
        """
        Adds the table wit the detected AMR mutations.
        :param section: Report section
        :return: None
        """
        data_mutations = pd.read_table(self._tool_inputs['TSV_point'][0].path)
        logger.info(f'{len(data_mutations)} mutations parsed')
        section.add_header('Detected AMR mutations', 3)
        section.add_table([[
            *[HtmlTableCell(f'{row[col]:.2f}' if isinstance(row[col], float) else row[col], color='green')
              for col in data_mutations.columns if col != 'PMID'],
            ResFinderReporter.__get_accession_cell(str(row['PMID']), 'green', is_pmid=True)
        ] for row in data_mutations.fillna('-').to_dict('records')], list(data_mutations.columns), [('class', 'data')])

        # Download link
        if len(data_mutations) > 0:
            relative_path = Path('resfinder4', self._tool_inputs['TSV_point'][0].path.name)
            section.add_file(self._tool_inputs['TSV_point'][0].path, relative_path)
            section.add_link_to_file('Download (TSV)', relative_path)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection(ResFinderReporter.TITLE, subtitle=self._input_informs['resfinder']['_name'])

        # Phenotype overviews
        for key in ('species', 'general'):
            if f'TSV_pheno_{key}' not in self._tool_inputs:
                continue
            self.__add_phenotype_table(section, key)
        section.add_warning_message(
            "The phenotype 'no resistance' should be interpreted with caution, as genes or mutations may be missing " 
            "from the database. In addition, these are WGS-based predictions that may not be reflected in the "
            "phenotype.")
        section.add_horizontal_line()

        # Genes & mutations
        self.__add_genes_table(section)
        if 'TSV_point' in self._tool_inputs:
            self.__add_mutations_table(section)
        section.add_horizontal_line()

        # Extra information
        self.__add_explanation_matches(section)
        section.add_paragraph(f"Database version: {self._input_informs['resfinder']['db_version_resfinder']}")
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]

    def __add_explanation_matches(self, section: HtmlReportSection) -> None:
        """
        Adds information about the different type of matches to the bottom of the report.
        :param section: Report section
        :return: None
        """
        section.add_header('Extra information', 3)
        section.add_paragraph('The following colors are used to denote the different type of hits:')
        section.add_table([
            [HtmlTableCell('', color='green'), 'Perfect match (100% over full length)'],
            [HtmlTableCell('', color='lightgreen'), 'Coverage 100%, identity <100%'],
            [HtmlTableCell('', color='grey'), 'Coverage <100%, identity <= 100%'],
            [HtmlTableCell('', color=None), 'No match found'],
        ], None, [('class', 'data')])
