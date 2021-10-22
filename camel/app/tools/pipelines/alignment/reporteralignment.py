from pathlib import Path
from typing import Union, List, Dict

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class ReporterAlignment(Tool):

    """
    Tool to create HTML reports for mapping
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the reporter
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Alignment: reporter', '0.1', camel)
        self._report_section = None
        self.__sub_folder = Path('alignment')

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - BAM, PDF_CG, PDF_MQC are required
        - For informs, alignment, picardmetrics, and samtoolsdepth are required.
        :return: None
        """
        for key in ['BAM', 'PDF_GC', 'PDF_MQC']:
            if key not in self._tool_inputs:
                raise InvalidInputSpecificationError(f'Required input {key} missing for alignment reporter')
            if len(self._tool_inputs[key]) != 1:
                raise InvalidInputSpecificationError(f'Illegal number of input files (max = 1) for input {key} '
                                                     f'provided for alignment reporter: {self._tool_inputs[key]}')

        for key in ['alignment', 'picardmetrics', 'samtoolsdepth']:
            if key not in self._input_informs:
                raise InvalidInputSpecificationError(f'Required informs missing for alignment reporter: {key}')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._report_section = HtmlReportSection('Alignment', subtitle=self._input_informs['alignment']['_name'])

        self._report_section.add_header('General', 3)
        self.__add_files()
        self._report_section.add_line_break()
        self.__add_alignment_metrics()

        self._report_section.add_header('Whole genome statistics', 3)
        self.__add_whole_genome_statistics()

        self._report_section.add_header('Segment statistics', 3)
        self.__add_segment_statistics()

        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section, False)]

    def __add_alignment_metrics(self) -> None:
        """
        Adds the general alignment metrics to the report.
        :return: None
        """
        table_data = [
            ['Total reads', self.__get_total_reads()],
            ['Aligned reads', self.__get_aligned_reads()],
            ['Reads aligned in pairs', self.__get_aligned_pairs()],
            ['Mean read length', self.__get_mean_read_length()],
            ['Median insert size', self.__get_median_insert_size()],
            ['Aligned bases', self.__get_aligned_bases()],
            ['Mismatch rate (per aligned base)', self.__get_mismatch_rate()],
            ['INDEL rate (per 100 aligned bases)', self.__get_indel_rate()],
            ['GC dropout (%)', self.__get_gc_dropout()]
        ]
        self._report_section.add_table(table_data, column_names=['Metric', 'Value'], table_attributes=[('class', 'data')])

    def __get_total_reads(self) -> str:
        """
        Returns the formatted string with the total number of reads.
        :return: Formatted string
        """
        return self.__reformat_inform(str(self._input_informs['picardmetrics']['AlignmentSummary_stats']['TOTAL_READS']))

    def __get_aligned_reads(self) -> str:
        """
        Returns the formatted string with the number of aligned reads.
        :return: Formatted string
        """
        reads = self._input_informs['picardmetrics']['AlignmentSummary_stats']['PF_READS_ALIGNED']
        pct = float(self._input_informs['picardmetrics']['AlignmentSummary_stats']['PCT_PF_READS_ALIGNED']) * 100
        return self.__reformat_inform(f'{reads} ({pct:.2f}%)')

    def __get_aligned_pairs(self) -> str:
        """
        Returns the formatted string with the number of aligned read pairs.
        :return: Formatted string
        """
        reads = self._input_informs['picardmetrics']['AlignmentSummary_stats']['READS_ALIGNED_IN_PAIRS']
        pct = float(self._input_informs['picardmetrics']['AlignmentSummary_stats']['PCT_READS_ALIGNED_IN_PAIRS']) * 100
        return self.__reformat_inform(f'{reads} ({pct:.2f}%)')

    def __get_mean_read_length(self) -> str:
        """
        Returns the formatted string with the mean read length.
        :return: Formatted string
        """
        return self.__reformat_inform(self._input_informs['picardmetrics']['AlignmentSummary_stats']['MEAN_READ_LENGTH'])

    def __get_median_insert_size(self) -> str:
        """
        Returns the formatted string with the median insert size.
        :return: Formatted string
        """
        return self.__reformat_inform(self._input_informs['picardmetrics']['InsertSize_stats']['MEDIAN_INSERT_SIZE'])

    def __get_aligned_bases(self) -> str:
        """
        Returns the formatted string with the number of aligned bases.
        :return: Formatted string
        """
        return self.__reformat_inform(self._input_informs['picardmetrics']['AlignmentSummary_stats']['PF_ALIGNED_BASES'])

    def __get_mismatch_rate(self) -> str:
        """
        Returns the formatted string with the mismatch rate.
        :return: Formatted string
        """
        mismatch = float(self._input_informs['picardmetrics']['AlignmentSummary_stats']['PF_MISMATCH_RATE']) * 100
        return f'{self.__reformat_inform(str(mismatch), 4)}%'

    def __get_indel_rate(self) -> str:
        """
        Returns the formatted string with the indel rate.
        :return: Formatted string
        """
        indel = float(self._input_informs['picardmetrics']['AlignmentSummary_stats']['PF_INDEL_RATE']) * 100
        return f'{self.__reformat_inform(str(indel), 4)}%'

    def __get_gc_dropout(self) -> str:
        """
        Returns the formatted string with the GC dropout.
        :return: Formatted string
        """
        at = float(self._input_informs['picardmetrics']['GCBias_stats']['AT_DROPOUT'])
        at_str = self.__reformat_inform(f'{at:.2f}')
        gc = float(self._input_informs['picardmetrics']['GCBias_stats']['GC_DROPOUT'])
        gc_str = self.__reformat_inform(f'{gc:.2f}')
        return f'GC {gc_str}, AT {at_str}'

    def __add_files(self) -> None:
        """
        Add the bam, bai and pdf files to the output report.
        :return: None
        """
        relative_path_bam = self.__sub_folder / 'mapped_reads.ba'
        self._report_section.add_file(self._tool_inputs['BAM'][0].path, f'{relative_path_bam}m')
        self._report_section.add_link_to_file('Alignment BAM file', f'{relative_path_bam}m')
        self._report_section.add_file(str(self._tool_inputs['BAM'][0].path).replace('.bam', '.bai'), f'{relative_path_bam}i')
        self._report_section.add_link_to_file('Alignment BAM index', f'{relative_path_bam}i')

        relative_path_pdf_gc = self.__sub_folder / 'gc_bias_curve.pdf'
        self._report_section.add_file(self._tool_inputs['PDF_GC'][0].path, str(relative_path_pdf_gc))
        self._report_section.add_link_to_file('GC bias curve (PDF)', str(relative_path_pdf_gc))

        relative_path_map_qc = self.__sub_folder / 'mapping_quality_distribution.pdf'
        self._report_section.add_file(self._tool_inputs['PDF_MQC'][0].path, str(relative_path_map_qc))
        self._report_section.add_link_to_file('Mapping quality distribution (PDF)', str(relative_path_map_qc))

    def __add_whole_genome_statistics(self) -> None:
        """
        Adds the table with the whole genome statistics to the report.
        :return: None
        """
        column_names = ['Sequence ID', 'Median cov.', 'Cov. MAD', 'Cov. IQR', 'Base cov. ratio (%)']
        depth_stats = self._input_informs['samtoolsdepth']
        table_data = [
            ['Whole genome',
             self.__reformat_inform(depth_stats['median_coverage'], 1),
             self.__reformat_inform(depth_stats['coverage_mad'], 1),
             self.__reformat_inform(depth_stats['coverage_iqr']),
             self.__reformat_inform(depth_stats['base_coverage'])]  # ratio already multipled with 100
        ]
        self._report_section.add_table(table_data, column_names=column_names, table_attributes=[('class', 'data')])

    def __add_segment_statistics(self) -> None:
        """
        Adds the statistics per identified segment to the output report. If no segments are specified, it only prints
        the sequence id. When segments are specified, it will give the results in the given segment order and with
        '-' when a segment was not identified in the analysis.
        :return: None
        """
        column_names = ['Sequence ID', 'Median cov.', 'Cov. MAD', 'Cov. IQR', 'Base cov. ratio (%)']
        stats = self._input_informs['samtoolsdepth']
        segments = self._parameters['genome_segments'].value.split(',') if 'genome_segments' in self._parameters else None
        table_data = []
        for segment, key in self.__get_segment_seqids(segments).items():
            if key is not None:
                table_data.append([key,
                                   self.__get_median_coverage_cell(key),
                                   self.__reformat_inform(stats['segment_coverage_mad'][key], 1),
                                   self.__reformat_inform(stats['segment_coverage_iqr'][key]),
                                   self.__get_base_cov_ratio_cell(key)
                                   ])
            else:
                table_data.append([f'({segment})', '-', '-', '-', '-'])
        self._report_section.add_table(table_data, column_names=column_names, table_attributes=[('class', 'data')])

    def __get_segment_seqids(self, segments: List[str]) -> Dict[str, str]:
        """
        Returns the sequence ids that have to be processed based on the segments. In case no segments were given, the
        sequence ids that are found are returned as both key and value. When segments are given, the key is the
        segment and the value is the sequence id corresponding to that segment.
        :param segments: List of segments
        :return: Dictionary with mapping between segment and sequence ids
        """
        if segments is None:
            if len(self._input_informs['samtoolsdepth']['segment_median_coverage']) != 1:
                return {k: k for k in self._input_informs['samtoolsdepth']['segment_median_coverage'].keys()}
        else:
            segment_order = {}
            for segment in segments:
                segment_order[segment] = None
                for key in self._input_informs['samtoolsdepth']['segment_median_coverage'].keys():
                    if key.endswith(f'-{segment}') or key.endswith(f'|{segment}'):
                        segment_order[segment] = key
            return segment_order

    def __get_median_coverage_cell(self, key: str) -> HtmlTableCell:
        """
        Returns a table cell with the median coverage. If the median coverage is below the specified
        threshold, the cell is colored orange.
        :param key: Sequence id to get the value for
        :return: Formatted table cell with median coverage
        """
        value = self._input_informs['samtoolsdepth']['segment_median_coverage'][key]
        if value < self._parameters['read_depth_cutoff'].value:
            return HtmlTableCell(self.__reformat_inform(value, 1), color='orange')
        else:
            return HtmlTableCell(self.__reformat_inform(value, 1))

    def __get_base_cov_ratio_cell(self, key: str) -> HtmlTableCell:
        """
        Returns a table cell with the base coverage ratio. If the median coverage is below the specified
        threshold, the cell is colored orange.
        :param key: Sequence id to get the value for
        :return: Formatted table cell with base coverage ratio
        """
        value = self._input_informs['samtoolsdepth']['segment_base_coverage'][key]
        if value < self._parameters['base_coverage_cutoff'].value:
            return HtmlTableCell(self.__reformat_inform(value, 2), color='orange')
        else:
            return HtmlTableCell(self.__reformat_inform(value, 2))

    @staticmethod
    def __reformat_inform(inp: Union[str, float, int], decimals: int = 2) -> str:
        """
        This function is used to reformat an inform value to a more readable format.
        This function also works when the percentage is omitted.
        E.g. 5241241 (10.02%) -> 5.241.241 (10.02%)
        :param inp: Input
        :param decimals: Decimal places to keep for floats
        :return: Reformatted inform
        """
        parts = str(inp).split(' ')
        if len(parts) == 1:
            try:
                return f'{int(parts[0]):,}'
            except ValueError:
                return f'{float(parts[0]):,.{decimals}f}'
        elif len(parts) == 2:
            try:
                return f'{int(parts[0]):,} {parts[1]}'
            except ValueError:
                return f'{float(parts[0]):,.{decimals}f} {parts[1]}'
        raise ValueError(f'Cannot parse: {inp}')
