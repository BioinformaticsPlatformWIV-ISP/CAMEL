import logging
import collections


STProfile = collections.namedtuple('STProfile', 'name alleles metadata')

from camel.app.tools.tool import Tool
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError


class SequenceTypeDetector(Tool):
    """
    Tool that manages MLST schemes. Also reports scheme metadata information in the informs.
    """

    SYMBOL_NO_ST = 'ND'

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: Camel instance
        :return: None
        """
        super(SequenceTypeDetector, self).__init__('Typing: Sequence Type Detector', '0.1', camel)
        self._wildcards = None
        self._symbol_allele_absent = None
        self._metadata_columns = []

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        profiles_file = self._tool_inputs['TSV'][0].path
        # allele_ids = {h.value.locus: h.value.allele_id for h in self._tool_inputs['VAL_Hits']}
        allele_ids = {}
        for h in self._tool_inputs['VAL_Hits']:
            if h.value.is_perfect_hit():
                allele_ids[h.value.locus] = h.value.allele_id
            else:
                allele_ids[h.value.locus] = '-'
        profiles = self.__parse_profiles(profiles_file, allele_ids.keys())
        sequence_type = self.__get_sequence_type(profiles, allele_ids)
        self._informs['sequence_type'] = sequence_type if sequence_type is not None else STProfile(
            SequenceTypeDetector.SYMBOL_NO_ST, None, [(k, '-') for k, _ in self._metadata_columns])
        logging.info("Detected sequence type: {}".format(self._informs['sequence_type'].name))

    def _check_input(self):
        """
        Checks whether the input is correct.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Sequence type profiles are required")
        if 'VAL_Hits' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No hits input found")
        super(SequenceTypeDetector, self)._check_input()

    def __parse_profiles(self, profiles_file, gene_names):
        """
        Parses the sequence type profiles.
        :param profiles_file: Profiles file
        :return: List of profiles
        """
        profiles = []
        with open(profiles_file) as handle:
            header = handle.readline().strip().split('\t')
            gene_indices = {gene_name: header.index(gene_name) for gene_name in gene_names}
            self._metadata_columns = [(p, header.index(p)) for p in header if p not in gene_names]
            for line in handle.readlines():
                parts = line.split('\t')
                parts[-1] = parts[-1].strip()
                alleles = {gene_name: parts[gene_indices[gene_name]] for gene_name in gene_names}
                metadata = [(name, parts[i]) for name, i in self._metadata_columns]
                profiles.append(STProfile(parts[0], alleles, metadata))
        return profiles

    def __alleles_match(self, detected_allele, profile_allele):
        """
        Checks whether two alleles match.
        :param detected_allele: Detected allele
        :param profile_allele: Allele from the ST profile
        :return: True if the alleles match
        """
        if profile_allele == self._symbol_allele_absent:
            return detected_allele == '-'
        elif profile_allele in self._wildcards:
            return True
        return detected_allele == profile_allele

    def __get_sequence_type(self, profiles, gene_alleles):
        """
        Returns the sequence type corresponding to the given gene alleles.
        :param profiles: Sequence type profiles
        :param gene_alleles: Alleles for the genes.
        :return: Sequence type and metadata
        """
        if 'allele_absent_symbol' in self._parameters:
            self._symbol_allele_absent = self._parameters['allele_absent_symbol'].value
        self._wildcards = self._parameters['allele_wildcard'].value.split(',')
        logging.debug("Wildcards: [{}]¸ Symbol allele absent: {}".format(
            ', '.join(self._wildcards), self._symbol_allele_absent))
        for profile in profiles:
            matched = True
            for gene_name, allele in gene_alleles.items():
                if not self.__alleles_match(allele, profile.alleles[gene_name]):
                    matched = False
                    break
            if matched:
                return profile
        return None
