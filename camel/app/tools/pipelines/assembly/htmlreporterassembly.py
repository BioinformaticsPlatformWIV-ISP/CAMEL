import os

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class HtmlReporterAssembly(Tool):
    """
    Tool to create HTML reports for the Assembly.
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('HTML Reporter', '0.1', camel)
        self.__subfolder = 'assembly'
        self._report_section = None

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self._report_section = HtmlReportSection('Assembly', subtitle=self._input_informs['spades']['_name'])
        self.__add_assembly_info()
        self.__add_assembly_download_link()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section, False)]

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if 'FASTA_Contig' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No assembly input found ('FASTA_Contig')")
        if 'SAMPLE_NAME' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No sample name input found ('SAMPLE_NAME')")
        if 'ASSEMBLER' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No assembler input found ('ASSEMBLER')")
        if 'quast' not in self._input_informs:
            raise InvalidInputSpecificationError("Quast informs are required ('quast')")
        if 'spades' not in self._input_informs:
            raise InvalidInputSpecificationError("SPAdes informs are required ('spades')")
        super()._check_input()

    def __add_assembly_info(self):
        """
        Adds the assembly info.
        :return: None
        """
        quast_informs = self._input_informs['quast']
        table_data = [
            ('Assembler:', self._tool_inputs['ASSEMBLER'][0].value),
            ('N50:', '{:,}'.format(int(quast_informs['contig']['N50']))),
            ('Number of contigs:', '{:,}'.format(int(quast_informs['contig']['# contigs (>= 1000 bp)']))),
            ('Total length:', '{:,}'.format(int(quast_informs['genome']['Total length'])))
        ]
        self._report_section.add_table(table_data, table_attributes=[('class', 'information')])

    def __add_assembly_download_link(self):
        """
        Adds a download link for the assembly.
        :return: None
        """
        relative_path = os.path.join(self.__subfolder, '{}_contigs.fasta'.format(
            FileSystemHelper.make_valid(self._tool_inputs['SAMPLE_NAME'][0].value)))
        self._report_section.add_file(self._tool_inputs['FASTA_Contig'][0].path, relative_path)
        self._report_section.add_link_to_file('Assembly (FASTA)', relative_path)
