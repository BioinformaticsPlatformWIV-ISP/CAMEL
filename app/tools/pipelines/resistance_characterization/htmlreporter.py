import os

from app.components.files.tsvexporter import TsvExporter
from app.components.html.tablecell import HtmlTableCell
from app.components.resistance.argannotparser import ArgannotParser
from app.components.resistance.cardparser import CardParser
from app.components.resistance.resfinderparser import ResFinderParser

from app.tools.export.htmlreporter import HtmlReporter


class HtmlReporterResistanceCharacterization(HtmlReporter):
    """
    Tool that creates HTML reports for the resistance characterization pipeline.
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: Camel instance
        :return: None
        """
        super(HtmlReporterResistanceCharacterization, self).__init__(camel)

    def _create_report(self):
        """
        Creates the HTML report for this pipeline.
        :return: None
        """
        self.__subfolder = os.path.join('resistance_characterization', self.__get_database_name().lower())
        self._report.add_header(self.__get_database_name(), 3)
        if len(self._tool_inputs['VAL_Hits']) == 0:
            self._report.add_paragraph('No hits found.')
        else:
            self.__add_gene_output_table()
        self._report.add_paragraph('Last updated: {}'.format(self._input_informs['database_info']['last_updated']))
        self._report.add_horizontal_line()
        self._tool_outputs['DIR'] = [self._tool_inputs['DIR'][0]]

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if 'database_info' not in self._input_informs:
            raise ValueError("No database info found")
        if 'VAL_Hits' not in self._tool_inputs:
            raise ValueError("No blast hits found")
        if 'TXT' not in self._tool_inputs:
            raise ValueError("No alignment inputs found")
        super(HtmlReporterResistanceCharacterization, self)._check_input()

    def __add_gene_output_table(self):
        """
        Adds the table with the detected genes.
        :return: None
        """
        header = ['Resistance gene', '%Identity', 'HSP/Gene length', 'Contig', 'Position in contig', 'Accession',
                  'Alignment']
        table_data = []
        for hit_input, txt_alignment in zip(self._tool_inputs['VAL_Hits'], self._tool_inputs['TXT']):
            hit = hit_input.value
            gene_name, annotation = self.__parse_header(hit.database_gene)
            row_color = hit.color
            table_data.append(
                [HtmlTableCell(gene_name, [('class', row_color)]),
                 HtmlTableCell(str(hit.percent_identity), [('class', row_color)]),
                 HtmlTableCell('{} / {}'.format(hit.alignment_length, hit.database_gene_length),
                               [('class', row_color)]),
                 HtmlTableCell(hit.query, [('class', row_color)]),
                 HtmlTableCell('{}..{}'.format(hit.query_start, hit.query_end),
                               [('class', row_color)]),
                 HtmlTableCell(annotation, [('class', row_color)],
                               'https://www.ncbi.nlm.nih.gov/nuccore/{}'.format(annotation)),
                 self.__save_alignment(txt_alignment, row_color)])
        self._report.add_table(table_data, header, [('class', 'data')])

        TsvExporter.export(table_data, header, os.path.join(self._folder, 'detected_genes.tmp'), [6])
        table_path = os.path.join(self.__subfolder, 'detected_genes.tsv')
        self._report.add_link_to_file('Download (TSV)',
                                      self._save_file(os.path.join(self._folder, 'detected_genes.tmp'), table_path))

    def __get_database_name(self):
        """
        Returns the name of the current database.
        :return: Name of the database
        """
        return self._input_informs['database_info']['name']

    def __parse_header(self, header):
        """
        Parses the header form the databases.
        :return: Gene name, annotation
        """
        database_name = self.__get_database_name()
        if database_name.lower() == 'argannot':
            header_data = ArgannotParser.parse_header(header)
            return header_data['gene_name'], header_data['accession']
        if database_name.lower() == 'card':
            header_data = CardParser.parse_header(header)
            return header_data['gene_name'], header_data['accession_gb']
        if database_name.lower() == 'resfinder':
            header_data = ResFinderParser.parse_header(header)
            return header_data['allele_name'], header_data['accession']
        else:
            raise ValueError("Unknown database: {}".format(database_name))

    def __save_alignment(self, alignment, cell_color):
        """
        Saves the alignment in the output folder.
        :param alignment: Alignment file
        :param cell_color: Color of the cell
        :return: Table cell containing a link to the alignment
        """
        path = os.path.join(self.__subfolder, 'alignments', os.path.basename(alignment.path))
        return HtmlTableCell('view', [('class', cell_color)], self._save_file(alignment.path, path))
