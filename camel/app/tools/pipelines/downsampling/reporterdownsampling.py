from typing import Any, Dict

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class ReporterDownsampling(Tool):
    """
    Creates output reports for the down sampling workflow.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Reporter: downsampling', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'stats' not in self._input_informs:
            raise InvalidInputSpecificationError("Statistics input is required")
        super()._check_input()

    def __add_table_params(self, section: HtmlReportSection, stats: Dict[str, Any]) -> None:
        """
        Adds the table with parameters.
        :param section: HTML report section
        :param stats: Dictionary with statistics
        :return: None
        """
        section.add_table([
            [f'Reference genome size:', f"{stats['size_ref_genome']:,} bases"],
            [f'Maximum coverage:', f"{stats['coverage_target']:.2f}"],
        ], table_attributes=[('class', 'information')])

    def __add_table_stats(self, section: HtmlReportSection, stats: Dict[str, Any], reads_are_paired: bool) -> None:
        """
        Adds the table with downsampling statistics.
        :param section: HTML report section
        :param stats: Dictionary with statistics
        :param reads_are_paired: True if reads are paired, False otherwise
        :return: None
        """
        text_downsample = 'No downsampling required' if (stats['downsample_factor'] is None) else \
            f"{stats['downsample_factor']:.2f}"
        section.add_table([
            [f"{int(stats['total_bases']):,}", str(int(stats['mean_read_length'])),
             f"{stats['coverage_estimated']:.2f}", text_downsample]],
            ['Total bases (fwd. + rev.)' if reads_are_paired else 'Total bases', 'Read length (avg.)',
             'Estimated coverage', 'Downsample factor'],
            [('class', 'data')]
        )

    def __add_table_downsampling(self, section: HtmlReportSection, stats: Dict[str, Any], reads_are_paired: bool) -> None:
        """
        Adds the table with the downsampling information.
        :param section: HTML report section
        :param stats: Dictionary with statistics
        :param reads_are_paired: True if reads are paired, False otherwise
        :return: None
        """
        import pprint
        pprint.pprint(stats)
        print(reads_are_paired)
        if reads_are_paired:
            header = ['Read pairs in', 'Read pairs out']
            table_data = [[f"{stats['nb_read_pairs_in']:,}", f"{self._input_informs['seqtk']['reads_count'] // 2:,}"]]
        else:
            header = ['Reads in', 'Reads out']
            table_data = [[f"{stats['nb_reads_in']:,}", f"{self._input_informs['seqtk']['reads_count']:,}"]]
        section.add_header_with_subtitle('Downsampling', 3, self._input_informs['seqtk']['_name'])
        section.add_table(table_data, header, [('class', 'data')])

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Create report section
        section = HtmlReportSection('Coverage check')
        reads_are_paired = 'is_paired' in self._parameters
        stats = self._input_informs['stats']

        # Add overview tables
        self.__add_table_params(section, stats)
        self.__add_table_stats(section, stats, reads_are_paired)

        # Add downsampling table if it was performed
        if ('seqtk' in self._input_informs) and (stats['downsample_factor'] is not None):
            self.__add_table_downsampling(section, stats, reads_are_paired)
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
