import os

from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class VariantCallingReporter(Tool):
    """
    This tool is used to create reports for the variant calling pipeline.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__("Variant Calling Reporter", "0.1", camel)
        self._section = HtmlReportSection("Variant calling")

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'VCF' not in self._tool_inputs:
            raise InvalidInputSpecificationError("VCF input is required")
        if 'VCF_filtered' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Filtered VCF input is required")
        if 'VAL_Sample' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Sample name is required")
        if 'reference' not in self._input_informs:
            raise InvalidInputSpecificationError("Reference information is required")
        if 'mapping' not in self._input_informs:
            raise InvalidInputSpecificationError("Mapping informs are required")
        super()._check_input()

    def __add_filtering_table(self):
        """
        Adds a table with the filtering information.
        :return: None
        """
        self._section.add_header('Filtering', 4)
        filtering_table_data = []
        for filter_name, filter_key in [
            ('Depth', 'filter_depth'),
            ('SNP quality', 'filter_snp_qual'),
            ('Mapping quality', 'filter_mapping_qual'),
            ('Distance', 'filter_distance'),
            ('Z-score', 'filter_zscore')
        ]:
            if filter_key in self._input_informs:
                stats = '{}/{}'.format(self._input_informs[filter_key]['variants_out'],
                                       self._input_informs[filter_key]['variants_in'])
            else:
                stats = 'NA'
            filtering_table_data.append((filter_name, stats,))
        self._section.add_table(filtering_table_data, ['Filter', 'Variants passed'], [('class', 'data')])

    def __add_output_files_table(self):
        """
        Adds a table containing the output files.
        :return: None
        """
        self._section.add_header('Output files', 4)
        vcf = self.__create_vcf_download_cell(self._tool_inputs['VCF'][0].path, False)
        vcf_filtered = self.__create_vcf_download_cell(self._tool_inputs['VCF_filtered'][0].path, True)
        table_data = [
            ['Unfiltered', self._input_informs['unfiltered']['total_variants'], vcf],
            ['Filtered', self._input_informs['filtered']['variants_out'], vcf_filtered]]
        self._section.add_table(table_data, ['', 'Number of variants', 'VCF file'], [('class', 'data')])

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__add_reference_link()
        self.__add_mapping_info()
        self.__add_filtering_table()
        self.__add_output_files_table()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]
        self.__set_informs()

    def __add_reference_link(self):
        """
        Adds a link to the reference genome that was used.;
        :return: None
        """
        info = self._input_informs['reference']
        reference_paragraph = HtmlElement('p', 'Reference: ')
        reference_paragraph.add_html_object(HtmlElement('a', info['name'], [('href', info['url'])]))
        self._section.add_html_object(reference_paragraph)

    def __add_mapping_info(self):
        """
        Adds the information about the read mapping to the report.
        :return: None
        """
        self._section.add_header('Read mapping', 4)
        table_data = [
            ['{}%'.format(self._input_informs['mapping']['stats_map_rate'])]
        ]
        self._section.add_table(table_data, ['Mapping rate'], [('class', 'data')])

    def __create_vcf_download_cell(self, path, filtered):
        """
        Creates a table cell that contains a download link for the VCF file.
        :param path: Path to the file
        :param filtered: True if this is the filtered VCF file
        :return: Table cell
        """
        if filtered:
            filename = f"variants-{self._tool_inputs['VAL_Sample'][0].value}-filtered.vcf"
        else:
            filename = f"variants-{self._tool_inputs['VAL_Sample'][0].value}.vcf"
        relative_path = os.path.join('variant_calling', filename)
        self._section.add_file(path, relative_path)
        return HtmlTableCell('Download', link=relative_path)

    def __set_informs(self):
        """
        Sets the informs to create the summary file.
        :return: None
        """
        self._informs.update({
            'depth_in': self._input_informs['filter_depth']['variants_in'],
            'depth_out': self._input_informs['filter_depth']['variants_out'],
            'dist_in': self._input_informs['filter_distance']['variants_in'],
            'dist_out': self._input_informs['filter_distance']['variants_out'],
            'mapping_qual_in': self._input_informs['filter_mapping_qual']['variants_in'],
            'mapping_qual_out': self._input_informs['filter_mapping_qual']['variants_out'],
            'snp_qual_in': self._input_informs['filter_snp_qual']['variants_in'],
            'snp_qual_out': self._input_informs['filter_snp_qual']['variants_out'],
            'zscore_in': self._input_informs['filter_zscore']['variants_in'],
            'zscore_out': self._input_informs['filter_zscore']['variants_out'],
            'mapping_rate': self._input_informs['mapping']['stats_map_rate']
        })
