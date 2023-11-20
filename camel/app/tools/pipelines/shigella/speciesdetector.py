import logging

import vcf

from camel.app.camel import Camel
from camel.app.tools.tool import Tool


class SpeciesDetector(Tool):
    """
    This tool is used to detect the species for the Shigella pipeline.
    Distinguishes E. coli and Shigella using
    - Presence / absence of ipaH gene
    - Presence / absence of TG indel in speG gene

    Input:
    - TSV: Tabular file containing the hits detected for the species identification gene database.

    Output:
    - INFORMS: 'detected_species' ('E. coli', 'Shigella' or 'NA')
    """

    SPEG_INDEL_POSITION = 2256002
    SPEG_START_POSITION = 2255997
    SPEG_END_POSITION = 2256035

    PROFILES = {
        'header': ['Species', '<i>ipaH</i>', '<i>speG</i> TG missing'],
        'data': [['<i>E. coli</i>', '-', '-'], ['<i>Shigella</i>', '+', '+']]
    }

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Shigella: species detector', '0.1', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._informs['ipaH_present'] = self.__check_presence_ipah()
        self._informs['speG_present'] = self.__check_presence_speg()
        if self._informs['speG_present'] is True:
            self._informs['speG_indel_present'] = self.__check_presence_speg_indel()
        else:
            self._informs['speG_indel_present'] = False
        self._informs['detected_species'] = self.__determine_species()

    def __check_presence_ipah(self) -> bool:
        """
        Checks if the ipaH gene is present (based on the gene detection hits).
        :return: True if present, False otherwise
        """
        hits = [io.value for io in self._tool_inputs['VAL_hits']]
        return 'ipaH' in [h.locus for h in hits]

    def __check_presence_speg(self) -> bool:
        """
        Checks if the speG gene is present (based on the samtools depth stats for the region).
        :return: True if present, False otherwise
        """
        speg_min_depth = int(self._parameters['threshold_speG_depth'].value)
        self._informs['speG_depth_threshold'] = speg_min_depth
        self._informs['speG_depth'] = self._input_informs['speG_depth']['median_depth']
        return self._input_informs['speG_depth']['median_depth'] > speg_min_depth

    def __check_presence_speg_indel(self) -> bool:
        """
        Checks if the speG TG indel is present by checking the VCF file at the given position.
        :return: True if present, False otherwise
        """
        with open(self._tool_inputs['VCF'][0].path) as handle:
            for variant in vcf.Reader(handle):
                if not (SpeciesDetector.SPEG_START_POSITION <= variant.POS <= SpeciesDetector.SPEG_END_POSITION):
                    continue
                logging.debug(f"Variant in speG: {variant}")
                if variant.POS == SpeciesDetector.SPEG_INDEL_POSITION and variant.ALT[0] == 'AGT':
                    return True
        return False

    def __determine_species(self) -> str:
        """
        Determines the species.
        :return: Species
        """
        if self._informs['ipaH_present'] is True:
            return 'Shigella'
        elif (self._informs['ipaH_present'] is False) and (self._informs['speG_indel_present'] is False):
            return 'E. coli'
        else:
            return 'NA'
