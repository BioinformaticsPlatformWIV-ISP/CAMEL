import os

from app.components.genedetection.genedetectionsrst2hit import GeneDetectionSRST2Hit
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue
from app.tools.tool import Tool


class SRST2HitExtractor(Tool):
    """
    This tool extracts hits from SRST2 output.

    INPUT:
        - TSV: SRST2 output file
        - mapping: Mapping of sequence ids to the original headers
        - db_info: Database information

    OUTPUT:
        - TSV: Tabular output file with the hits
        - VAL_Hits: List os SRST2 Hit objects
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Gene Detection: SRST2 Hit Extractor', '0.1', camel)

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        hits = []
        if len(self._tool_inputs['TSV']) > 0:
            tsv_input = self._tool_inputs['TSV'][0].path
            with open(tsv_input) as handle:
                for line in handle.readlines()[1:]:
                    hits.append(GeneDetectionSRST2Hit.create_from_srst2_output_line(
                        line, self._input_informs['mapping'], self._input_informs['db_info']['metadata']))
        self._tool_outputs['VAL_Hits'] = sorted([ToolIOValue(h) for h in hits], key=lambda v: v.value.locus)
        output_path = os.path.join(self._folder, self._parameters['output_filename'].value)
        self.__create_output_file(hits, output_path)
        self._tool_outputs['TSV'] = [ToolIOFile(output_path)]

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("TSV input is required")
        if 'mapping' not in self._input_informs:
            raise InvalidInputSpecificationError("Sequence id mapping is required")
        if 'db_info' not in self._input_informs:
            raise InvalidInputSpecificationError("Database information is required")
        super()._check_input()

    @staticmethod
    def __create_output_file(hits, path):
        """
        Creates the tabular output file.
        :param hits: Detected hits
        :param path: Output path
        :return: None
        """
        with open(path, 'w') as handle:
            if len(hits) > 0:
                handle.write('\t'.join(hits[0].get_table_column_names()))
                handle.write('\n')
                for hit in hits:
                    handle.write(hit.to_table_row())
                    handle.write('\n')
