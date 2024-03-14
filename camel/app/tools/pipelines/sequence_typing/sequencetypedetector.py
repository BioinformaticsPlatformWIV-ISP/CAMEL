import logging
import pandas as pd
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
from pathlib import Path

from camel.app.components.sequencetyping.sequencetypinghitbase import SequenceTypingHitBase
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.tools.tool import Tool
from camel.app.io.tooliofile import ToolIOFile


@dataclass(frozen=True, unsafe_hash=True)
class STProfile:
    name: str
    metadata: List[Tuple[str, str]] = field(hash=False)
    alleles: Optional[Dict[str, str]] = field(hash=False)


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
        super().__init__('Typing: Sequence Type Detector', '0.1', camel)
        self._wildcards = None
        self._symbol_allele_absent = None
        self._metadata_columns = []

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse input data
        hit_by_locus = {hit.locus: hit for hit in [
            io.value for io in self._tool_inputs['hits_nucl'] + self._tool_inputs['hits_pept']]}
        profiles = self.__parse_profiles(list(hit_by_locus.keys()))

        # Determine best profile
        nb_matches_by_profile = self.__get_nb_matches_by_profile(profiles, hit_by_locus)
        best_profile = sorted(nb_matches_by_profile.items(), key=lambda x: -x[-1])[0][0]
        percent_matching = 100 * nb_matches_by_profile[best_profile] / len(best_profile.alleles)
        is_detected = percent_matching >= int(self._parameters['min_percent_detected'].value)

        # Write all profiles to a tsv file
        if ('write_tsv' in self._parameters.keys()) and self._parameters['write_tsv'].value:
            num_alleles = len(list(nb_matches_by_profile.keys())[0].alleles)
            all_matches = pd.DataFrame([(stprofile.name, hits / num_alleles) for stprofile, hits in
                                        nb_matches_by_profile.items()], columns=["ST", "proportion_match"])
            all_matches.to_csv(self._folder / Path(self._parameters['output_filename'].value), sep="\t", index=False)
            self._tool_outputs['TSV_all_matches'] = [ToolIOFile(self._folder / Path(self._parameters['output_filename'].value))]

        # Save output data
        self._informs.update({
            'is_detected': is_detected,
            'nb_detected': nb_matches_by_profile[best_profile],
            'nb_loci': len(best_profile.alleles),
            'percent_detected': percent_matching,
            'symbol': best_profile.name if is_detected else SequenceTypeDetector.SYMBOL_NO_ST,
            'metadata': [(k, v if is_detected else '-') for k, v in best_profile.metadata]
        })

    def _check_input(self) -> None:
        """
        Checks whether the input is correct.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Sequence type profiles ('TSV') are required.")
        if len(self._tool_inputs['hits_nucl']) + len(self._tool_inputs['hits_pept']) == 0:
            raise InvalidInputSpecificationError("Typing hits are required.")
        super(SequenceTypeDetector, self)._check_input()

    def __parse_profiles(self, gene_names: List[str]) -> List[STProfile]:
        """
        Parses the sequence type profiles.
        :param gene_names: Name of the genes
        :return: List of profiles
        """
        profiles = []
        with open(self._tool_inputs['TSV'][0].path) as handle:
            header = handle.readline().strip().split('\t')
            gene_indices = {gene_name: header.index(gene_name) for gene_name in gene_names}
            self._metadata_columns = [(p, header.index(p)) for p in header if p not in gene_names]
            for line in handle.readlines():
                parts = line.split('\t')
                parts[-1] = parts[-1].strip()
                alleles = {gene_name: parts[gene_indices[gene_name]] for gene_name in gene_names}
                metadata = [(name, parts[i]) for name, i in self._metadata_columns]
                profiles.append(STProfile(name=parts[0], alleles=alleles, metadata=metadata))
        return profiles

    def __alleles_match(self, detected_hit: SequenceTypingHitBase, profile_allele: str) -> bool:
        """
        Checks whether two alleles match.
        :param detected_hit: Detected allele
        :param profile_allele: Allele from the ST profile
        :return: True if the alleles match
        """
        if profile_allele == self._symbol_allele_absent:
            return detected_hit.allele_id == '-'
        elif profile_allele in self._wildcards:
            return True
        return (detected_hit.is_perfect_hit()) and (detected_hit.allele_id == profile_allele)

    def __get_nb_matches_by_profile(
            self, profiles: List[STProfile], hit_by_locus: Dict[str, SequenceTypingHitBase]) -> Dict:
        """
        Returns the number of matches per profile.
        :param profiles: Sequence type profiles
        :param hit_by_locus: Alleles for the genes.
        :return: Nb. of matches by profile
        """
        if 'allele_absent_symbol' in self._parameters:
            self._symbol_allele_absent = self._parameters['allele_absent_symbol'].value
        self._wildcards = self._parameters['allele_wildcards'].value.split(',')
        logging.debug("Wildcards: [{}]¸ Symbol allele absent: {}".format(
            ', '.join(self._wildcards), self._symbol_allele_absent))

        nb_matches_by_profile = {}
        for profile in profiles:
            nb_matching = sum([
                self.__alleles_match(hit, profile.alleles[gene_name]) for gene_name, hit in hit_by_locus.items()])
            nb_matches_by_profile[profile] = nb_matching
        return nb_matches_by_profile

# # TODO: remove
# if __name__ == '__main__':
#     from camel.app.camel import Camel
#     from camel.app.snakemake.snakemakeutils import SnakemakeUtils
#     camel = Camel.get_instance()
#
#     ## Enterocolitica
#     hits_nucl = Path('/home/fistrijthaegen/PycharmProjects/camel_yersinia/test/enterocolitica/typing/cgmlst/DNA/hits.io')
#     hits_pept = Path('/home/fistrijthaegen/PycharmProjects/camel_yersinia/test/enterocolitica/typing/cgmlst/peptide/hits.io')
#     tsv = Path('/home/fistrijthaegen/PycharmProjects/camel_yersinia/test/enterocolitica/typing/cgmlst/tsv-profiles.io')
#
#     # Pseudotuberculosis
#     # hits_nucl = Path('/home/fistrijthaegen/PycharmProjects/camel_yersinia/test/pseudotuberculosis/typing/cgmlst/DNA/hits.io')
#     # hits_pept = Path('/home/fistrijthaegen/PycharmProjects/camel_yersinia/test/pseudotuberculosis/typing/cgmlst/peptide/hits.io')
#     # tsv = Path('/home/fistrijthaegen/PycharmProjects/camel_yersinia/test/pseudotuberculosis/typing/cgmlst/tsv-profiles.io')
#
#     # Neisseria meningitidis
#     # hits_nucl = Path('/home/fistrijthaegen/PycharmProjects/camel_yersinia/test/neisseria/typing/cgmlst/DNA/hits.io')
#     # hits_pept = Path('/home/fistrijthaegen/PycharmProjects/camel_yersinia/test/neisseria/typing/cgmlst/peptide/hits.io')
#     # tsv = Path('/home/fistrijthaegen/PycharmProjects/camel_yersinia/test/neisseria/typing/cgmlst/tsv-profiles.io')
#
#     detector = SequenceTypeDetector(camel)
#     SnakemakeUtils.add_pickle_inputs(detector, {'hits_nucl':hits_nucl, 'hits_pept': hits_pept,'TSV':tsv})
#     detector.update_parameters(write_tsv='True', output_filename='profile_matches.tsv')
#     detector.run()
#     SnakemakeUtils.dump_tool_output(detector, 'TSV_matches', Path('tsv_profile_matches.io'))