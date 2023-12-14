from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool
import json
from pathlib import Path


class SPIFinderReporter(Tool):
    """
    Parses SPIFinder csv output reports and creates an html report.
    """

    TITLE = 'SPIFinder'

    def __init__(self, camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('SPIFinder Reporter', '0.1', camel)
        self._section = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        if 'spifinder_fastq' in self._input_informs:
            self._section = HtmlReportSection(SPIFinderReporter.TITLE,
                                              subtitle=self._input_informs['spifinder_fastq']['_name'])
        else:
            self._section = HtmlReportSection(SPIFinderReporter.TITLE,
                                              subtitle=self._input_informs['spifinder_fasta']['_name'])
        self._section.add_header('HITS - KMA on raw reads (FASTQ)', 3)
        if 'JSONFASTQ' in self._tool_inputs:
            self.__add_hits_results(self._tool_inputs['JSONFASTQ'][0].path, 'fastq')
        else:
            self._section.add_paragraph('SPIFinder raw reads results not available in FASTA-input mode')
        self._section.add_horizontal_line()
        self._section.add_header('HITS - BLAST on the assembly (FASTA)', 3)
        self.__add_hits_results(self._tool_inputs['JSONFASTA'][0].path, 'fasta')
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]
        self.__add_output_table_link()
        relative_path = Path('spifinder', 'summary_out.tsv')
        self._section.add_file(self._tool_inputs['TSV_output'][0].path, relative_path)
        relative_path_doc = Path('spifinder', 'spifinder_function_category.tsv')
        self._section.add_link_to_file("Category function definition table (TSV)", relative_path_doc)
        self._section.add_file(self._tool_inputs['TSV_doc'][0].path, relative_path_doc)
        self.__add_database_information()

    def __add_hits_results(self, res_file_path: Path, mode: str) -> None:
        """
        Adds the table with the antibiotic sensitivity
        :param res_file_path: path to the results file
        :param mode: either fasta or fastq
        :return: None
        """

        header = ['SPI', 'identity', 'HSP/locus length', 'Contig', 'Positions in contig', 'Accession', 'Insertion site',
                  'Category function']
        data = []

        if mode == 'fasta':
            column_to_keep = list(range(8))
            col_ncbi = 5
        else:  # mode == 'fastq':
            column_to_keep = [0, 1, 2, 5, 6, 7]
            col_ncbi = 3

        with open(res_file_path) as json_file:
            handle = json.load(json_file)
            spi = handle['spifinder']["results"]['Salmonella Pathogenicity Islands']['SPI']
            if spi == "No hit found":
                self._section.add_paragraph('No hits found.')
            else:
                for hits in spi.keys():
                    if spi[hits]['identity'] == 100:
                        color = 'green'
                    elif spi[hits]['identity'] >= 95:
                        color = 'lightgreen'
                    else:
                        color = 'yellow'
                    row = []
                    input_cols = [spi[hits]['SPI'], f"{spi[hits]['identity']:.2f}",
                                  f"{spi[hits]['HSP_length']}/{spi[hits]['template_length']}" if spi[hits].get('HSP_length') else spi[hits]['coverage'],
                                  spi[hits]['contig_name'] if spi[hits].get('contig_name') else '',
                                  spi[hits]['positions_in_contig'] if spi[hits].get('contig_name') else '',
                                  spi[hits]['accession'], spi[hits]['insertion_site'], spi[hits]['category_function']]
                    input_cols = [input_cols[index] for index in column_to_keep]
                    for i in range(len(input_cols)):
                        if i == col_ncbi:  # if it's the accession number then add link
                            link = f'https://www.ncbi.nlm.nih.gov/nuccore/{input_cols[i]}'
                            row.append(HtmlTableCell(input_cols[i], color, link=link))
                        else:
                            row.append(HtmlTableCell(input_cols[i], color))
                    data.append(row)
                self._section.add_table(data, [header[index] for index in column_to_keep], [('class', 'data')])

    def __add_output_table_link(self) -> None:
        """
        add output table link
        @return: None
        """
        relative_path = Path('spifinder', 'summary_out.tsv')
        self._section.add_link_to_file("Download (TSV)", relative_path)

    def __add_database_information(self) -> None:
        """
        Adds the database information to the report.
        :return: None
        """
        # spifinder_fasta is always executed and therefore this db_update_date is used
        self._section.add_paragraph('Last updated: {}'.format(self._input_informs['spifinder_fasta'].get(
            'last_update_date', '{LAST_UPDATE_DATE}')))
