from typing import Optional, List, Dict

import yaml

from camel.app.camel import Camel
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class SCCmecTyping(Tool):
    """
    This tool is used to run the SCCmecFinder tool.
    """
    COMPLEXES = [
        {'key': 'ccr_genes_complexes', 'name': 'ccr gene complex'},
        {'key': 'mec_genes_complexes', 'name': '<i>mec</i> gene complex'},
        {'key': 'SCC_mec_types', 'name': 'SCC<i>mec</i> type'}
    ]

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('SCCmec typing', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided tool input is valid.
        :return: None
        """
        if 'YML' not in self._tool_inputs:
            raise ValueError("Profiles input is required (YML)")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse profiles
        with open(self._tool_inputs['YML'][0].path) as handle:
            profile_data = yaml.load(handle, Loader=yaml.SafeLoader)

        # Get the detected genes
        detected_genes = [hit.locus.split(':')[0] for hit in [io.value for io in self._tool_inputs['VAL_HITS']]]

        # Determine the complexes
        self._informs['complexes'] = []
        for complex_ in SCCmecTyping.COMPLEXES:
            matching_complex = SCCmecTyping.__get_matching_complex(detected_genes, profile_data[complex_['key']])
            self._informs['complexes'].append({'value': matching_complex, **complex_})

    @staticmethod
    def __get_matching_complex(detected_genes: List[str], genes_by_complex: Dict[str, List[str]]) -> \
            Optional[str]:
        """
        Returns the matching complex (if there is one).
        :param genes_by_complex: Genes by complex
        :return: Complex (or None if there is none found)
        """
        for complex_, genes in genes_by_complex.items():
            if all(g in detected_genes for g in genes):
                return complex_
        logger.debug("No complex found")
