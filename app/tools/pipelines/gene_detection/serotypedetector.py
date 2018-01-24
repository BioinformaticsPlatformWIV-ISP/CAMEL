from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliovalue import ToolIOValue
from app.tools.tool import Tool


class SerotypeDetector(Tool):
    """
    This tool detects the serotype based on gene hits.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Serotype Detector', '0.1', camel)

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        for char in ('H', 'O'):
            if 'VAL_Hits_{}'.format(char) not in self._tool_inputs:
                raise InvalidInputSpecificationError("{}-type hits are required".format(char))
            if 'metadata_{}'.format(char) not in self._input_informs:
                raise InvalidInputSpecificationError("{}-type metadata informs are required".format(char))
        super()._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        h_type = self.__get_h_type()
        o_type = self.__get_o_type()
        self._tool_outputs['VAL_H_type'] = [ToolIOValue(h_type)]
        self._tool_outputs['VAL_O_type'] = [ToolIOValue(o_type)]
        self._tool_outputs['VAL_serotype'] = [ToolIOValue('{}:{}'.format(h_type, o_type))]

    def __get_h_type(self):
        """
        Returns the H type of the input.
        :return: H type.
        """
        metadata = self._input_informs['metadata_H']['metadata']
        hits_h = [v.value.subject for v in self._tool_inputs['VAL_Hits_H']]
        if len(hits_h) == 0:
            return 'H-ambiguous'
        elif len(hits_h) == 1:
            return metadata[hits_h[0]].get('predicted_serotype', '-')
        else:
            hits_sorted = {'fliC': [h for h in hits_h if 'fliC' in h],
                           'non-fliC': [h for h in hits_h if 'fliC' not in h]}
            if len(hits_sorted['non-fliC']) == 1:
                return metadata[hits_sorted['non-fliC'][0]].get('predicted_serotype', '-')
        return 'H-ambiguous'

    def __get_o_type(self):
        """
        Returns the O type.
        :return: O type.
        """
        metadata = self._input_informs['metadata_O']['metadata']
        hits_o = [v.value.subject for v in self._tool_inputs['VAL_Hits_O']]
        if len(hits_o) == 0:
            return 'O-ambiguous'
        elif len(hits_o) == 1:
            return metadata[hits_o[0]].get('predicted_serotype', '-')
        else:
            hits_sorted = {'xy': [h for h in hits_o if any([gene in h for gene in ('wzx', 'wzy')])],
                           'mt': [h for h in hits_o if any([gene in h for gene in ('wzm', 'wzt')])]}
            if len(hits_sorted['xy']) > 0 and len(hits_sorted['mt']) > 0:
                return 'O-ambiguous'
            else:
                all_hits = hits_sorted['xy'] + hits_sorted['mt']
                o_type = metadata[all_hits[0]].get('predicted_serotype', '-')
                for hit in all_hits:
                    if metadata[hit].get('predicted_serotype', '-') != o_type:
                        return "O-ambiguous"
                return o_type
