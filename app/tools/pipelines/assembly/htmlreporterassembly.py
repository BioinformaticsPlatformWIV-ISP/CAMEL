import os

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.tools.export.htmlreporter import HtmlReporter


class HtmlReporterAssembly(HtmlReporter):
    """
    Tool to create HTML reports for the Assembly.
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: CAMEL instance
        :return: None
        """
        super(HtmlReporterAssembly, self).__init__(camel)
        self.__subfolder = 'assembly'

    def _create_report(self):
        """
        Creates the HTML report.
        :return: None
        """
        self._report.add_header('Assembly', 2)
        if 'velvet_optimiser' in self._input_informs:
            self.__add_velvetoptimiser_info()
        else:
            self.__add_spades_info()
        self.__add_assembly_download_link()
        self._report.add_horizontal_line()

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if 'FASTA_Contig' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No assembly input found")
        if 'SAMPLE_NAME' not in self._tool_inputs:
            raise ValueError("No sample name input found")
        super(HtmlReporterAssembly, self)._check_input()

    def __add_velvetoptimiser_info(self):
        """
        Adds a table with the VelvetOptimiser info.
        :return: None
        """
        informs = self._input_informs['velvet_optimiser']
        table_data = [
            ['Assembler:', 'VelvetOptimiser'],
            ['Kmer used:', informs['kmer']],
            ['N50:', informs['n50']]
        ]
        self._report.add_table(table_data, table_attributes=[('class', 'information')])

    def __add_spades_info(self):
        """
        Adds the SPAdes info.
        :return: None
        """
        table_data = [
            [['Assembler:', 'SPAdes'],
             ['Kmer user:', '/'],
             ['N50:', '/']]
        ]
        self._report.add_table(table_data, table_attributes=[('class', 'information')])

    def __add_assembly_download_link(self):
        """
        Adds a download link for the assembly.
        :return: None
        """
        assembly_path = os.path.join(self.__subfolder, '{}_contigs.fasta'.format(
            self._tool_inputs['SAMPLE_NAME'][0].value))
        self._save_file(self._tool_inputs['FASTA_Contig'][0].path, assembly_path)
        self._report.add_link_to_file('Assembly (FASTA)', assembly_path)
