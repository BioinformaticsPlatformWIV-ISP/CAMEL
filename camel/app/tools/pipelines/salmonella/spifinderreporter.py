import json
from pathlib import Path

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class SPIFinderReporter(Tool):
    """
    Parses SPIFinder csv output reports and creates a html report.
    """

    TITLE = 'SPIFinder'

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('SPIFinder Reporter', '0.1', camel)
        self._section = None
        self._fastq_results_present = True if 'spifinder_fastq' in self._input_informs else False

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        super(SPIFinderReporter, self)._check_input()
        if self._fastq_results_present and 'JSON_FASTQ' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Fastq analysis results were found in the input informs; "
                                                 "JSON_FASTQ is required as input for this tool.")
        if 'JSON_FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("JSON_FASTA is missing but always required as input for this tool.")
        if 'TSV_output' not in self._tool_inputs:
            raise InvalidInputSpecificationError("TSV_output is missing but always required as input for this tool.")
        if 'TSV_documentation' not in self._tool_inputs:
            raise InvalidInputSpecificationError("TSV_documentation is missing but always required as input for "
                                                 "this tool.")

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._section = HtmlReportSection(SPIFinderReporter.TITLE,
                                          subtitle=self._input_informs['spifinder_fasta']['_name'])
        # Add Fastq results 'section'
        self._section.add_header('HITS - KMA on raw reads (FASTQ)', 3)
        if self._fastq_results_present:
            self.__add_hits_results(self._tool_inputs['JSON_FASTQ'][0].path, 'fastq')
        else:
            self._section.add_paragraph('SPIFinder raw reads results not available in FASTA-input mode')
        self._section.add_horizontal_line()
        # Add mandatory Fasta results 'section'
        self._section.add_header('HITS - BLAST on the assembly (FASTA)', 3)
        self.__add_hits_results(self._tool_inputs['JSON_FASTA'][0].path, 'fasta')

        self.__add_file_output()
        self.__add_database_information()

        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_hits_results(self, res_file_path: Path, mode: str) -> None:
        """
        Adds the table with the antibiotic sensitivity.
        :param res_file_path: path to the results file
        :param mode: either fasta or fastq
        :return: None
        """
        table_header = ['SPI', 'identity', 'HSP/locus length', 'Contig', 'Positions in contig', 'Accession',
                        'Insertion site', 'Category function']

        # Fasta gives more meaningful output fields than Fastq; the not meaningful ones are filtered out
        if mode == 'fasta':
            column_to_keep = list(range(len(table_header)))
            col_ncbi = table_header.index('Accession')
        elif mode == 'fastq':
            column_to_keep = [0, 1, 2, 5, 6, 7]
            table_header_subset = [table_header[index] for index in column_to_keep]
            col_ncbi = table_header_subset.index('Accession')
        else:
            raise ValueError(f"This function's parameter 'mode' must be either fastq or fasta, current value is {mode}")

        with open(res_file_path) as handle:
            json_file = json.load(handle)
        spi = json_file['spifinder']["results"]['Salmonella Pathogenicity Islands']['SPI']
        if spi == "No hit found":
            self._section.add_paragraph('No hits found.')
            return
        table_data = []
        for hit in spi.keys():
            color = self.___assign_hit_color(spi[hit]['identity'])
            hit_data = [spi[hit]['SPI'],  # SPI
                          f"{spi[hit]['identity']:.2f}",  # identity
                          f"{spi[hit]['HSP_length']}/{spi[hit]['template_length']}" if spi[hit].get('HSP_length')
                          else spi[hit]['coverage'],  # coverage
                          spi[hit].get('contig_name', ''),  # contig
                          spi[hit].get('positions_in_contig', ''),  # positions in contig
                          spi[hit]['accession'],  # accession
                          spi[hit]['insertion_site'],  # insertion site
                          spi[hit]['category_function']]  # category function

            # Keep columns based on input mode
            hit_data = [hit_data[index] for index in column_to_keep]

            # Add all values to the row
            row = [HtmlTableCell(value, color) for value in hit_data]

            # Add a href (blue clickable url redirect) to the ncbi accession column by overwriting it
            row[col_ncbi] = HtmlTableCell(hit_data[col_ncbi], color,
                                          link=f'https://www.ncbi.nlm.nih.gov/nuccore/{hit_data[col_ncbi]}')

            table_data.append(row)

        self._section.add_table(table_data, [table_header[index] for index in column_to_keep], [('class', 'data')])

    def ___assign_hit_color(self, percentage_identity: int) -> str:
        """
        Assigns a color according to percentage thresholds.
        :param percentage_identity: percentage identity for a hit.
        :return: color as str
        """
        if percentage_identity == 100:
            return 'green'
        elif percentage_identity >= 95:
            return 'lightgreen'
        else:
            return 'yellow'

    def __add_database_information(self) -> None:
        """
        Adds the database information to the report.
        :return: None
        """
        # spifinder_fasta is always executed and therefore this db_update_date is used
        self._section.add_paragraph('Last updated: {}'.format(self._input_informs['spifinder_fasta'].get(
            'last_update_date', '{LAST_UPDATE_DATE}')))

    def __add_file_output(self) -> None:
        """
        Add the output tsv files to the output report section.
        :return: None
        """
        relative_path = Path('spifinder', 'summary_out.tsv')
        self._section.add_link_to_file("Download (TSV)", relative_path)
        self._section.add_file(self._tool_inputs['TSV_output'][0].path, relative_path)
        relative_path_doc = Path('spifinder', 'spifinder_function_category.tsv')
        self._section.add_link_to_file("Category function definition table (TSV)", relative_path_doc)
        self._section.add_file(self._tool_inputs['TSV_documentation'][0].path, relative_path_doc)
