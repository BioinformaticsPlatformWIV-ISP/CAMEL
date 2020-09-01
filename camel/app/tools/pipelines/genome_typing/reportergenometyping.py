from pathlib import Path

from camel.app.camel import Camel
from camel.app.components.html.htmlexpandabletable import HtmlExpandableTable
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool
from typing import List


class ReporterGenomeTyping(Tool):
    """
    Tool to create HTML reports for genome typing.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize this tool.
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('Genome Typing: reporter', '0.1', camel)
        self._report_section = None
        self.__sub_folder = Path('genometyping')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._report_section = HtmlReportSection('De novo genome segment typing')
        self.__add_genome_typing_section()
        self.__add_file_output()
        self.__add_subtype_section()
        self.__add_segment_sections()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section, False)]

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'genometyping' not in self._input_informs:
            raise InvalidInputSpecificationError("No genome typing info found")
        super()._check_input()

    def __add_genome_typing_section(self) -> None:
        """
        Adds the segment typing section.
        :return: None
        """
        table_data = [
            ['Genome segments expected', ','.join(self._input_informs['genometyping']['expected_segments'])],
            ['Genome segment(s) recovered', self.__get_segment_string(self._input_informs['genometyping']['segment_coverage']['segment_covered'])],
            ['Genome segment(s) missing', self.__get_segment_string(self._input_informs['genometyping']['segment_coverage']['segment_missing'])],
            ['Number of reads for subtyping', self.__reformat_inform(str(self._input_informs['seqtksubsample']['reads_count']))]
        ]
        self._report_section.add_table(table_data, table_attributes=[('class', 'data')])

    @staticmethod
    def __get_segment_string(segments: List[str]):
        return ','.join(segments) if len(segments) != 0 else '-'

    def __add_file_output(self) -> None:
        """
        Saves the reference genome sequences and blast output to the output directory and adds links.
        :return: None
        """
        blast_filename = 'genometyping_blast.tsv'
        relative_path = self.__sub_folder / blast_filename
        self._report_section.add_file(self._tool_inputs['TSV'][0].path, str(relative_path))
        self._report_section.add_link_to_file('Detailed typing BLASTn results (TSV)', str(relative_path))

        ref_filename = 'genome_reference_segments.fasta'
        relative_path = self.__sub_folder / ref_filename
        self._report_section.add_file(self._tool_inputs['FASTA'][0].path, str(relative_path))
        self._report_section.add_link_to_file('Refence genome segment sequences (FASTA)', str(relative_path))

        self._report_section.add_line_break()
        self._report_section.add_horizontal_line()

    def __add_segment_sections(self) -> None:
        """
        Adds the sections per segment.
        :return: None
        """
        self._report_section.add_header('Typing results per segment', 2)
        self._report_section.add_text('Only the top five hits are shown, click to expand to see the top twenty hits')
        total_reads_count = self._input_informs['seqtksubsample']['reads_count']
        for i, segment in enumerate(self._input_informs['genometyping']['expected_segments']):

            self._report_section.add_header(f'Results on segment {segment}', 3)
            if segment in self._input_informs['genometyping']['segment_coverage']['segment_covered']:
                best_candidate_table = [['Best reference:', self._input_informs['genometyping']['segment_informs'][segment]['refseqid']]]
                if len(self._input_informs['genometyping']['segment_informs'][segment]['candidates']) > 1:
                    best_candidate_table.append(['Candidates:', ','.join(self._input_informs['genometyping']['segment_informs'][segment]['candidates'])])
            else:
                best_candidate_table = [['Best reference:', 'Segment not found in reads!']]
            self._report_section.add_table(best_candidate_table, table_attributes=[('class', 'information')])
            self._report_section.add_line_break()

            candidate_table = []
            for cnt in self._input_informs['genometyping']['segment_informs'][segment]['counts'][:20]:
                candidate_table.append([cnt[0], f'{self.__reformat_inform(str(cnt[1]))} ({cnt[1]/total_reads_count*100:.2f}%)'])
            self._report_section.add_html_object(HtmlExpandableTable(candidate_table, ['Segment name', 'Reads count (percentage)']))

            self._report_section.add_line_break()
            if i + 1 < len(self._input_informs['genometyping']['expected_segments']):
                self._report_section.add_horizontal_line()

    def __add_subtype_section(self) -> None:
        """
        Adds a section for the Influenza A subtyping results.
        :return: None
        """
        if 'hana_subtyping' in self._input_informs['genometyping']:
            informs = self._input_informs['genometyping']['hana_subtyping']
            if informs['failure_message']:
                self._report_section.add_warning_message(informs['failure_message'])
            self._report_section.add_header('Influenza A subtype detection', 2)
            table = [['Subtype', informs['subtype']],
                     ['Hemagglutinin (HA) subtype', informs['ha']],
                     ['Neuraminidase (NA) subtype', informs['na']]]
            self._report_section.add_table(table, table_attributes=[('class', 'data')])
            self._report_section.add_horizontal_line()

    @staticmethod
    def __reformat_inform(input_str: str) -> str:
        """
        This function is used to reformat an inform value to a more readable format.
        This function also works when the percentage is omitted.
        E.g. 5241241 (10.02%) -> 5.241.241 (10.02%)
        :param input_str: Input string
        :return: Reformatted inform
        """
        parts = input_str.split(' ')
        if len(parts) == 1:
            return f'{int(parts[0]):,}'
        elif len(parts) == 2:
            return f'{int(parts[0]):,} {parts[1]}'
        raise ValueError(f"Cannot parse: {input_str}")
