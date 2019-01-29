import collections
from typing import List

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class SerotypeDetectorEcoli(Tool):
    """
    This tool detects the E. coli serotype based on gene hits.
    """

    SerotypeHit = collections.namedtuple('SerotypeHit', 'gene type_')

    def __init__(self, camel: Camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Serotype Detector', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        for char in ('H', 'O'):
            if 'TSV_{}'.format(char) not in self._tool_inputs:
                raise InvalidInputSpecificationError("{}-type hits are required".format(char))
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        h_type = self.__get_h_type()
        o_type = self.__get_o_type()
        self._tool_outputs['VAL_H_type'] = [ToolIOValue(h_type)]
        self._tool_outputs['VAL_O_type'] = [ToolIOValue(o_type)]
        self._tool_outputs['VAL_serotype'] = [ToolIOValue('{}:{}'.format(o_type, h_type))]

    @staticmethod
    def __parse_hits(path) -> List[SerotypeHit]:
        """
        Parses the hits from a tabular gene detection output file.
        :param path: Path
        :return: List of hits
        """
        detected_genes = []
        with open(path) as handle:
            for line in handle.readlines()[1:]:
                parts = line.split('\t')
                # noinspection PyCallByClass
                detected_genes.append(SerotypeDetectorEcoli.SerotypeHit(parts[0], parts[-2]))
        return detected_genes

    def __get_h_type(self) -> str:
        """
        Returns the H type based on the gene detection hits.
        :return: H type.
        """
        detected_genes = SerotypeDetectorEcoli.__parse_hits(self._tool_inputs['TSV_H'][0].path)
        if len(detected_genes) == 0:
            return 'H-ambiguous'
        elif len(detected_genes) == 1:
            return detected_genes[0].type_
        else:
            hits_sorted = {'fliC': [h for h in detected_genes if 'fliC' in h.gene],
                           'non-fliC': [h for h in detected_genes if 'fliC' not in h.gene]}
            if len(hits_sorted['non-fliC']) == 1:
                return hits_sorted['non-fliC'][0].type_
        return 'H-ambiguous'

    def __get_o_type(self) -> str:
        """
        Returns the O type.
        :return: O type.
        """
        detected_genes = SerotypeDetectorEcoli.__parse_hits(self._tool_inputs['TSV_O'][0].path)
        if len(detected_genes) == 0:
            return 'O-ambiguous'
        elif len(detected_genes) == 1:
            return detected_genes[0].type_
        else:
            hits_sorted = {'xy': [h for h in detected_genes if any([g in h.gene for g in ('wzx', 'wzy')])],
                           'mt': [h for h in detected_genes if any([g in h.gene for g in ('wzm', 'wzt')])]}
            if len(hits_sorted['xy']) > 0 and len(hits_sorted['mt']) > 0:
                return 'O-ambiguous'
            else:
                all_hits = hits_sorted['xy'] + hits_sorted['mt']
                o_type = all_hits[0].type_
                for hit in all_hits:
                    if hit.type_ != o_type:
                        return "O-ambiguous"
                return o_type
