import json
from pathlib import Path
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool
from camel.app.camel import Camel


class SerovarSalmonellaReporter(Tool):
    """
    Parses Sistr tab output reports.
    """

    TITLE = 'Serotyping'

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Serovar Salmonella Reporter', '0.1', camel)
        self._section = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._section = HtmlReportSection(SerovarSalmonellaReporter.TITLE,
                                          subtitle=self._input_informs['serotyping_sistr']['_name'])
        with open(self._tool_inputs['TSV_SISTR'][0].path) as handle:
            data_sistr = json.load(handle)[0]
            if data_sistr['qc_status'] == 'PASS':
                self._section.add_header('H1 (fliC)-type', 4)
                self.__add_table_serotype_sistr_h_type('h1')
                self._section.add_header('H2 (fljB)-type', 4)
                self.__add_table_serotype_sistr_h_type('h2')
                self._section.add_header('O-type', 4)
                self.__add_table_serotype_sistr_h_type('o')
                self._section.add_header('Conclusion', 4)
                self.__add_conclusion()
            else:
                self._section.add_paragraph('SISTR did not pass QC: the isolate might be contaminated, '
                                            'or not contain enough sequencing reads, or ...')
        input_folder = self._tool_inputs['DIR_sistr'][0].path
        self.__add_database_information(input_folder)

        self._section.add_line_break()
        self._section.add_horizontal_line()
        self._section.add_header(self._input_informs['serotyping_seqsero2']['_name'], 2)
        self._section.add_header('SeqSero2 serotyping - raw read allele mode', 4)
        if 'TXTSeqSero2allele' in self._tool_inputs:
            self.__add_table_serotype_seqsero(self._tool_inputs['TXTSeqSero2allele'][0].path)
        else:
            self._section.add_paragraph('SeqSero2 serotyping (raw read allele mode) not available in FASTA-input mode')
        self._section.add_header('SeqSero2 serotyping - raw read kmer mode', 4)
        if 'TXTSeqSero2kmerread' in self._tool_inputs:
            self.__add_table_serotype_seqsero(self._tool_inputs['TXTSeqSero2kmerread'][0].path)
        else:
            self._section.add_paragraph('SeqSero2 serotyping (assembly kmer mode) not available in FASTA-input mode')
        self._section.add_header('SeqSero2 serotyping - assembly kmer mode', 4)
        self.__add_table_serotype_seqsero(self._tool_inputs['TXTSeqSero2kmer'][0].path)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]
        self.__add_output_table_link()
        # add file to the report
        relative_path = Path('serotyping', 'summary_out.tsv')
        self._section.add_file(self._tool_inputs['TSV_output'][0].path, relative_path)
        # add last update of the seqsero2 db
        input_folder2 = self._tool_inputs['DIR_seqsero2'][0].path
        self.__add_database_information(input_folder2)

    def __add_table_serotype_sistr_h_type(self, h: str) -> None:
        """
        Generates and adds the table of h1 type for sistr tool.
        :return: None
        """
        with open(self._tool_inputs['TSV_SISTR'][0].path) as json_file:
            handle = json.load(json_file)[0]
            table_data = []
            header = ['Locus', '% Identity', 'HSP/Locus length', 'Contig', 'Position in contig', 'Predicted serotype']
        if h == 'h1':
            prediction = [handle['h1_flic_prediction']]
            locus_full = ['fliC']
            locus_short = ['h1']
        elif h == 'h2':
            prediction = [handle['h2_fljb_prediction']]
            locus_full = ['fljB']
            locus_short = ['h2']
        else:  # h == 'o':
            prediction = [handle['serogroup_prediction']['wzx_prediction'],
                          handle['serogroup_prediction']['wzy_prediction']]
            locus_full = ['wzx', 'wzy']
            locus_short = ['serogroup']
        matcher = []
        for pred in range(len(prediction)):
            if str(prediction[pred]['is_missing']) == 'False':
                matcher.append(True)
                # col flic_prediction[pred]
                if prediction[pred]['top_result']['pident'] == 100.0 and \
                        prediction[pred]['top_result']['length'] == \
                        prediction[pred]['top_result']['qlen']:
                    color = 'green'
                elif prediction[pred]['top_result']['length'] == prediction[pred]['top_result']['qlen']:
                    color = 'lightgreen'
                else:
                    color = 'grey'

                row = [HtmlTableCell(x, color) for x in
                       [locus_full[pred], format(prediction[pred]['top_result']['pident'], '.2f'),
                        '/'.join([str(prediction[pred]['top_result']['length']),
                                  str(prediction[pred]['top_result']['qlen'])]),
                        prediction[pred]['top_result']['stitle'],
                        '...'.join([str(prediction[pred]['top_result']['sstart']),
                                    str(prediction[pred]['top_result']['send'])]),
                        prediction[pred][locus_short[0]]
                        ]
                       ]
                table_data.append(row)

            else:
                matcher.append(False)
                if h == 'o':
                    row = [HtmlTableCell(locus_full[pred], 'red'),
                           HtmlElement('td', 'No match found', attributes=[('colspan', 5), ('class', 'red')])]
                    table_data.append(row)
        if any(matcher):
            self._section.add_table(table_data, header, [('class', 'data')])
        if h == 'o':
            self._section.add_paragraph(' '.join(['Predicted O antigen based on H antigens and serogroup:',
                                                  handle['o_antigen']]))
        if not any(matcher):
            if h == 'h1' or h == 'h2':
                self._section.add_paragraph(f'{locus_full[0]}: No match found')
            else:
                self._section.add_paragraph(f'{locus_full[0]} and {locus_full[1]}: No match found')

    def __add_conclusion(self) -> None:
        """
        Generates and adds the conclusion sentence for sistr tool.
        :return: None
        """
        with open(self._tool_inputs['TSV_SISTR'][0].path) as json_file:
            handle = json.load(json_file)[0]
            self._section.add_paragraph(' '.join(['Predicted antigenic profile (O:H1:H2):',
                                                  ':'.join([str(handle['o_antigen']),
                                                            str(handle['h1']),
                                                            str(handle['h2'])])]))
            self._section.add_paragraph(' '.join(['Predicted serotype:',
                                                  handle['serovar']]))

    def __add_output_table_link(self) -> None:
        """
        Add the hyperlink for the tsv table for those assays (serotyping)
        :return: None
        """
        relative_path = Path('serotyping', 'summary_out.tsv')
        self._section.add_link_to_file("Download (TSV)", relative_path)

    def __add_table_serotype_seqsero(self, input_file_path: Path) -> None:
        """
        Generates and adds the table for seqsero2 tool.
        :return: None
        """
        with open(input_file_path) as txt_file:
            lines = txt_file.readlines()
            resultsdict = {}
            for line in lines:
                parts = line.split('\t')
                resultsdict[parts[0]] = parts[1]
            table_data = []
            header = ['O-antigen', 'H1-antigen (fliC)', 'H2-antigen (fljB)', 'Antigenic formula', 'Serotype']
            row = [resultsdict['O antigen prediction:'], resultsdict['H1 antigen prediction(fliC):'],
                   resultsdict['H2 antigen prediction(fljB):'], resultsdict['Predicted antigenic profile:'],
                   resultsdict['Predicted serotype:']]
            table_data.append(row)
            self._section.add_table(table_data, header, [('class', 'data')])

    def __add_database_information(self, input_folder: Path) -> None:
        """
        Adds the date of latest database update.
        :param input_folder: Input database directory
        :return: None
        """
        path_metadata = input_folder / 'db_update_info.json'
        if not path_metadata.is_file():
            raise FileNotFoundError(f'Database metadata not found: {path_metadata}')
        with path_metadata.open() as handle:
            metadata = json.load(handle)
            last_update = metadata['last_update_date']
        self._section.add_paragraph(f'Last updated: {last_update}')
