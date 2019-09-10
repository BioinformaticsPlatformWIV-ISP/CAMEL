from camel.app.camel import Camel
from camel.app.tools.tool import Tool


class SubspeciesDetector(Tool):
    """
    This tool is used to detect the sub-species for the Shigella pipeline.

    Input:
    - TSV: Tabular file containing the hits detected for the species identification gene database.

    Output:
    - VAL_subspecies: Detected sub-species
    """

    PROFILES = {
        'header': ['Subspecies'] + ['<i>{}</i>'.format(x) for x in [
            'rfc', 'wbgZ', 'wzx/wzy boydii', 'wzx/wzy dysenteriae']],
        'data': [
            ['<i>Flexneri</i>', '+', '-', '-', '-'],
            ['<i>Sonnei</i>', '-', '+', '-', '-'],
            ['<i>Boydii</i>', '-', '-', '+', '-'],
            ['<i>Dysenteriae</i>', '-', '-', '-', '+']
        ]
    }

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Shigella: subspecies detector', '0.1', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        hits = [io.value for io in self._tool_inputs['VAL_hits']]

        # ipaH and speG are used for the species, not for the subspecies
        detected_loci = [h.locus for h in hits if h.locus not in ('ipaH', 'speG')]
        self._informs['detected_loci_subspecies'] = ', '.join(sorted(detected_loci))

        # Determine the subspecies
        if all(['dysenteriae' in n for n in detected_loci]) and len(detected_loci) >= 1:
            self._informs['detected_subspecies'] = 'dysenteriae'
        elif all(['boydii' in n for n in detected_loci]) and len(detected_loci) >= 1:
            self._informs['detected_subspecies'] = 'boydii'
        elif len(detected_loci) == 1 and detected_loci[0] == 'rfc':
            self._informs['detected_subspecies'] = 'flexneri'
        elif len(detected_loci) == 1 and detected_loci[0] == 'wbgZ':
            self._informs['detected_subspecies'] = 'sonnei'
        else:
            self._informs['detected_subspecies'] = 'NA'
