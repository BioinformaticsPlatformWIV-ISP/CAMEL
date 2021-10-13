from pathlib import Path
from typing import List

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class ReporterConsensus(Tool):
    """
    Tool to create HTML reports for consensus sequence extraction.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize this tool.
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('Consensus: reporter', '0.1', camel)
        self._report_section = None
        self.__sub_folder = Path('conseq')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._report_section = HtmlReportSection('Consensus sequence extraction')
        self.__add_file_output()
        self.__add_segment_sections()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section, False)]

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'stats' not in self._input_informs:
            raise InvalidInputSpecificationError("No genome typing info found")
        super()._check_input()

    @staticmethod
    def __get_segment_string(segments: List[str]):
        return ','.join(segments) if len(segments) != 0 else '-'

    def __add_file_output(self) -> None:
        """
        Saves the reference genome sequences and blast output to the output directory and adds links.
        :return: None
        """
        sequence_variants = 'sequence_variants.vcf'
        relative_path = self.__sub_folder / sequence_variants
        self._report_section.add_file(self._tool_inputs['VCF'][0].path, str(relative_path))
        self._report_section.add_link_to_file('Sequence(s) variants (VCF)', str(relative_path))

        bam_filename = 'bam_processed.bam'
        relative_path = self.__sub_folder / bam_filename
        self._report_section.add_file(self._tool_inputs['BAM'][0].path, str(relative_path))
        self._report_section.add_link_to_file('Variant caller pre-processed input file (BAM)', str(relative_path))

        bai_filename = 'bam_processed.bai'
        relative_path = self.__sub_folder / bai_filename
        self._report_section.add_file(str(self._tool_inputs['BAM'][0].path).replace('.bam', '.bai'), str(relative_path))
        self._report_section.add_link_to_file('Variant caller input BAM index file (BAI)', str(relative_path))

        consensus_filename = 'consensus.fasta'
        relative_path = self.__sub_folder / consensus_filename
        self._report_section.add_file(self._tool_inputs['FASTA'][0].path, str(relative_path))
        self._report_section.add_link_to_file('Consensus sequence(s) (FASTA)', str(relative_path))

        blastn_filename = 'blastn_genome_reference_segments.tsv'
        relative_path = self.__sub_folder / blastn_filename
        self._report_section.add_file(self._tool_inputs['TSV'][0].path, str(relative_path))
        self._report_section.add_link_to_file('Blastn against consensus sequence (TSV)', str(relative_path))

        self._report_section.add_line_break()
        self._report_section.add_horizontal_line()

    def __add_segment_sections(self) -> None:
        """
        Adds the sections per segment.
        :return: None
        """
        self._report_section.add_header('Consensus sequence(s) statistics', 2)
        column_names = ['start', 'end', 'length', 'SNP count', 'N count', 'base cov. ratio (%)', 'median cov.', 'RefSeq length', 'insertions (count/bases)', 'deletions (count/bases)']
        column_ids_basic = ['start', 'end', 'length', 'snp_count', 'n_count', 'coverage', 'median_coverage', 'query_length']
        column_ids_special = ['insertion_count', 'inserted_base_count', 'deletion_count', 'deleted_base_count']
        informs = self._input_informs['stats']['target_consensus_seq_inform']
        for i, identifier in enumerate(informs.keys()):
            self._report_section.add_header(f'Results on {identifier}', 3)
            table_data = []
            for j, hit in enumerate(informs[identifier]):
                informs[identifier][j]['median_coverage'] = self._input_informs['samtools']['segment_median_coverage'][identifier]
                table_data_hit = [informs[identifier][j][x] for x in column_ids_basic]
                table_data_hit.append(f"{informs[identifier][j][column_ids_special[0]]}/{informs[identifier][j][column_ids_special[1]]}")
                table_data_hit.append(f"{informs[identifier][j][column_ids_special[2]]}/{informs[identifier][j][column_ids_special[3]]}")
                table_data.append(table_data_hit)
            self._report_section.add_table(table_data, column_names=column_names, table_attributes=[('class', 'data')])
