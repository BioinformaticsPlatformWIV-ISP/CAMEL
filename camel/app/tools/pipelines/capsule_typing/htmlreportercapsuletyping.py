import datetime

from camel.app.components.html.tablecell import HtmlTableCell
from camel.app.tools.export.htmlreporter import HtmlReporter
from camel.app.tools.pipelines.sequence_typing.besthitextractor import BestHitExtractor


class HtmlReporterCapsuleTyping(HtmlReporter):
    """
    Tool to create HTML reports for the capsule typing.
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
        :param camel: CAMEL instance
        :return: None
        """
        super(HtmlReporterCapsuleTyping, self).__init__(camel)
        self.__subfolder = 'capsule_typing'

    def _create_report(self):
        """
        Creates the HTML report.
        :return: None
        """
        self._report.add_header('Serogroup determination', 2)
        self._create_output_table()
        self._add_last_updated()
        self._report.add_horizontal_line()

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        super(HtmlReporterCapsuleTyping, self)._check_input()
        if 'typing' not in self._input_informs:
            raise ValueError("No capsule typing info found")

    @staticmethod
    def __get_detected_genes(serogroup_info):
        """
        Returns the detected.
        :param serogroup_info: Serogroup information ({'gene': {info}})
        :return: Gene names, number of detected genes
        """
        gene_names = []
        nb_of_detected_genes = 0
        for gene_name in sorted(serogroup_info.keys()):
            gene_names.append(gene_name)
            if serogroup_info[gene_name]['hit_type'] == BestHitExtractor.PERFECT_HIT:
                nb_of_detected_genes += 1
        return gene_names, nb_of_detected_genes

    def _create_output_table(self):
        """
        Creates the output table.
        :return: None
        """
        typing_informs = self._input_informs['typing']
        header = ['Serogroup', 'Genes detected'] + [''] * 7
        table_data = []
        for serogroup in sorted(typing_informs.keys()):
            gene_names, nb_of_detected_genes = HtmlReporterCapsuleTyping.__get_detected_genes(typing_informs[serogroup])
            row = [serogroup, '{0:.0f}%'.format(100 * float(nb_of_detected_genes) / len(gene_names))]
            for i in range(0, 7):
                try:
                    gene_name = gene_names[i]
                    color = HtmlReporterCapsuleTyping.COLOR_CODES[typing_informs[serogroup][gene_name]['hit_type']]
                    row.append(HtmlTableCell(gene_name, attributes=[('class', color)]))
                except IndexError:
                    row.append('-')
            table_data.append(row)
            table_data.sort(key=lambda x: float(x[1].replace('%', '')), reverse=True)
        self._report.add_table(table_data, header, [('class', 'data')])

    def _add_last_updated(self):
        """
        Adds the date when the database was last updated to the report.
        :return: None
        """
        try:
            date = datetime.date(*[int(x) for x in self._input_informs['last_updated'].split('-')]).strftime('%d-%m-%Y')
        except ValueError:
            date = 'NA'
        self._report.add_paragraph('Last updated: {}'.format(date))
