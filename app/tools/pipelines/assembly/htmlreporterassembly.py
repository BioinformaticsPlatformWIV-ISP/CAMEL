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
        self.__add_assembly_info(self._tool_inputs['ASSEMBLER'][0].value)
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
            raise InvalidInputSpecificationError("No sample name input found")
        if 'ASSEMBLER' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No assembler input found")
        if 'quast' not in self._input_informs:
            raise InvalidInputSpecificationError("Quast informs are required")
        super(HtmlReporterAssembly, self)._check_input()

    def __add_assembly_info(self, assembler):
        """
        Adds the assembly info.
        :param assembler: Name of the assembler that was used.
        :return: None
        """
        quast_informs = self._input_informs['quast']
        table_data = [
            ('Assembler', assembler),
            ('N50', quast_informs['contig'].get('N50', '-')),
            ('Number of contigs', quast_informs['contig'].get('# contigs (>= 0 bp)', '-')),
            ('Number of contigs (>1000bp)', quast_informs['contig'].get('# contigs (>= 1000 bp)', '-')),
            ('Total length', quast_informs['genome'].get('Total length'))
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
