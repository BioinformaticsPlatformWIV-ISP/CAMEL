import datetime
import os

from app.components.filesystemhelper import FileSystemHelper
from app.components.html.tablecell import HtmlTableCell
from app.tools.export.htmlreporter import HtmlReporter
from app.tools.pipelines.sequence_typing.besthitextractor import BestHitExtractor


class HtmlReporterSequenceTyping(HtmlReporter):
    """
    Tool that creates HTML reports for the sequence typing pipeline.
    """

    COLOR_CODES = {
        BestHitExtractor.NO_HIT: 'red',
        BestHitExtractor.MULTIPLE_HITS: 'yellow',
        BestHitExtractor.PERFECT_HIT: 'green',
        BestHitExtractor.IMPERFECT_IDENTITY_HIT: 'lightgreen',
        BestHitExtractor.IMPERFECT_LENGTH_HIT: 'grey'
    }

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: Camel instance
        :return: None
        """
        super(HtmlReporterSequenceTyping, self).__init__(camel)
        self.__subfolder = None

    def _create_report(self):
        """
        Creates the HTML report for this pipeline.
        :return: None
        """
        self.__subfolder = os.path.join('sequence_typing', self.__get_locus_set_directory_name().lower())
        self._report.add_header(self.__scheme_title(), 3)
        if self._input_informs['locus_info']['has_profile_definitions']:
            self.__add_sequence_type_table()
        self.__add_output_table()
        self.__add_output_download_link()
        self.__add_database_info()
        self.__add_linked_data()
        self._report.add_horizontal_line()

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if 'DIR_LS' not in self._tool_inputs:
            raise ValueError("No locus set directory input found")
        if 'TSV' not in self._tool_inputs:
            raise ValueError("No allele detection table input found")
        if 'locus_info' not in self._input_informs:
            raise ValueError("No locus set manager info found")
        if 'hits' not in self._input_informs:
            raise ValueError("No hit detection info found")
        if 'metadata_alleles' not in self._input_informs:
            raise ValueError("No allele metadata info found")
        if 'metadata_sequence_type' not in self._input_informs:
            raise ValueError("No sequence type metadata found")
        if 'SAMPLE_NAME' not in self._tool_inputs:
            raise ValueError("No sample name input found")
        super(HtmlReporterSequenceTyping, self)._check_input()

    def __get_locus_set_directory_name(self):
        """
        Returns the directory of the current locus set.
        :return: Directory of the locus set
        """
        return os.path.splitext(self._tool_inputs['DIR_LS'][0].basename)[0]

    def __scheme_title(self):
        """
        Returns the title for the scheme.
        :return: Scheme title
        """
        return self._input_informs['locus_info']['scheme_metadata']['title']

    def __add_sequence_type_table(self):
        """
        Adds the sequence type metadata.
        :return: None
        """
        metadata = self._input_informs['metadata_sequence_type']['metadata']
        sequence_type_cell = HtmlTableCell(metadata.values()[0], [('class', self.__get_sequence_type_color())])
        table_data = [[sequence_type_cell] + metadata.values()[1:]]
        header = [key.replace('_', ' ') for key in metadata.keys()]
        self._report.add_table(table_data, header, [('class', 'data')])

    def __get_sequence_type_color(self):
        """
        Returns the color for the sequence type cell.
        :return: None
        """
        if self._input_informs['sequence_type_detection']['sequence_type'] == 'ND':
            return 'red'
        else:
            return 'green'

    def __add_output_table(self):
        """
        Adds the output table to the report.
        :return: None
        """
        table_data = []
        header = ['Locus', 'Allele', '% Identity', 'HSP length / Locus length', 'Type', 'Alignment']
        with open(self._tool_inputs['TSV'][0].path) as input_handle:
            for line in input_handle.readlines()[1:]:
                locus, allele_id, identity, length, type_, hit = line.split('\t')
                table_row = [locus, self.__get_allele_id_cell(allele_id, locus, hit.strip()), identity, length, type_,
                             self.__save_alignment(locus)]
                table_data.append(table_row)
        self._report.add_table(table_data, header, table_attributes=[('class', 'data')])

    def __get_allele_id_cell(self, id_, locus, type_):
        """
        Returns the a HtmlTableCell for the given allele id.
        :param id_: Allele id
        :param locus: Locus
        :param type_: Hit type
        :return: Table cell
        """
        url = self._input_informs['hits'][locus]['url']
        color = HtmlReporterSequenceTyping.COLOR_CODES[type_]
        return HtmlTableCell(id_, [('class', color)], url)

    def __save_alignment(self, locus):
        """
        Saves the alignment of the given locus in the output folder.
        :param locus: Name of the locus
        :return: Table cell containing a link to the alignment
        """
        for alignment in self._tool_inputs.get('TXT', []):
            if alignment.basename == '{}.txt'.format(locus):
                path = os.path.join(self.__subfolder, 'alignments', alignment.basename)
                return HtmlTableCell('view', link=self._save_file(alignment.path, path))
        return '-'

    def __add_output_download_link(self):
        """
        Adds a download link for the output table.
        :return: None
        """
        table_path = os.path.join(self.__subfolder, 'alleles-{}-{}.tsv'.format(
            self.__get_locus_set_directory_name().lower(),
            FileSystemHelper.make_valid(self._tool_inputs['SAMPLE_NAME'][0].value)))
        self._save_file(self._tool_inputs['TSV'][0].path, table_path)
        self._report.add_link_to_file('Download (TSV)', table_path)

    def __add_database_info(self):
        """
        Adds database information.
        :return: None
        """
        informs = self._input_informs['locus_info']['scheme_metadata']
        date = datetime.date(*[int(x) for x in informs['last_updated'].split('-')])
        self._report.add_paragraph('Last update: {}'.format(date.strftime('%d-%m-%Y')))
        self._report.add_line_break()

    def __add_linked_data(self):
        """
        Adds the linked metadata to the report.
        :return: None
        """
        linked_data = self._input_informs['metadata_alleles'].copy()
        if len(linked_data) == 0:
            return
        table_data = []
        for allele_name in sorted(linked_data.keys()):
            label, value = linked_data[allele_name]
            table_data.append(['{}:'.format(label), value])
        self._report.add_table(table_data, table_attributes=[('class', 'information')])
