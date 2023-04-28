import json

import os
import vcf
from Bio import SeqIO

from camel.app.camel import Camel
from camel.app.components.csq.mutations.frameshiftmutation import FrameshiftMutation
from camel.app.components.csq.mutations.stopmutation import StopMutation
from camel.app.tools.tool import Tool


COORD_MINUS10_TA_BOX = [-13, -4]
COORD_MINUS35_BOX = [-30, -38]
START_GTR_OPERON = 65


class FlexneriTypeDetector(Tool):
    """
    This tool is used to detect the Shigella Flexneri type.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance.
        :return: None
        """
        super().__init__('Shigella: flexneri type detector', '0.1', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        locus_by_allele = self.__parse_fasta_input()
        self._informs['mutations'] = {}
        self._informs['loci'] = {}
        for allele, locus in locus_by_allele.items():
            print(locus)
            if locus not in self._informs['loci']:
                self._informs['loci'][locus] = {'detected': False, 'mutations': {'stop': [], 'frameshift': []}}
            key = f'VAL_mut_{allele}'
            if key not in self._tool_inputs:
                continue
            self._informs['loci'][locus]['detected'] = True
            self._informs['loci'][locus]['mutations']['stop'] = [
                io.value.long_notation for io in self._tool_inputs[key] if isinstance(io.value, StopMutation)]
            self._informs['loci'][locus]['mutations']['frameshift'] = [
                io.value.long_notation for io in self._tool_inputs[key] if isinstance(io.value, FrameshiftMutation)]
            self._informs['loci'][locus]['expressed'] = all([len(ms) == 0 for x, ms in self._informs['loci'][locus][
                'mutations'].items()])
            self._informs['loci'][locus]['VCF'] = self._tool_inputs[f'VCF_csq_{allele}'][0].path
        self.__parse_gtr_promotor_mutations(self._tool_inputs['VCF'][0].path)
        # Determine the type based on the profiles
        self.__parse_profiles(self._tool_inputs['TSV'][0].path)
        self._informs['detected_type'] = self.__get_flexneri_type()

    def __parse_fasta_input(self):
        """
        Parses the input FASTA file.
        :return:
        """
        locus_by_allele = {}
        with open(os.path.join(self._tool_inputs['DIR_FASTA'][0].path, 'shigella_flexneri_type.fasta')) as handle:
            for seq in SeqIO.parse(handle, 'fasta'):
                metadata = json.loads(' '.join(seq.description.split(' ')[1:]))
                locus_by_allele[metadata['allele']] = metadata['locus']
        return locus_by_allele

    def __parse_profiles(self, profiles_path: str) -> None:
        """
        Parses the profiles for the flexneri type.
        :param profiles_path: Profiles path
        :return: None
        """
        self._informs['profiles'] = {}
        with open(profiles_path) as handle:
            self._informs['all_loci'] = handle.readline().strip().split('\t')[1:]
            for line in handle.readlines():
                parts = line.strip().split('\t')
                self._informs['profiles'][parts[0]] = {
                    locus: True if parts[i+1] == '+' else False for i, locus in enumerate(self._informs['all_loci'])
                }

    def __get_flexneri_type(self) -> str:
        """
        Returns the Flexneri type, if none matches 'NA' is returned.
        :return: Flexneri type
        """
        for profile, state_by_locus in self._informs['profiles'].items():
            if all([self._informs['loci'][l].get('expressed', False) == s for l, s in state_by_locus.items()]):
                return profile
        return 'NA'

    def __parse_gtr_promotor_mutations(self, vcf_file: str) -> None:
        """
        Parses the mutations in the gtr promotor region.
        :param vcf_file: Input VCF file
        :return: None
        """
        with open(vcf_file) as handle:
            vcf_records = list(vcf.VCFReader(handle))

        self._informs['promotor_variants'] = {'-35_box': [], '-10_TA_box': [], 'other': []}
        self._informs['wt_gtr_promotor'] = True
        for variant in vcf_records:
            relative_position = variant.POS - START_GTR_OPERON + 1
            if relative_position > 0:
                continue
            elif COORD_MINUS35_BOX[0] <= relative_position <= COORD_MINUS35_BOX[1]:
                self._informs['promotor_variants']['-35_box'].append([relative_position, variant])
                self._informs['wt_gtr_promotor'] = False
            elif COORD_MINUS10_TA_BOX[0] <= relative_position <= COORD_MINUS10_TA_BOX[1]:
                self._informs['promotor_variants']['-10_TA_box'].append([relative_position, variant])
                self._informs['wt_gtr_promotor'] = False
            else:
                self._informs['promotor_variants']['other'].append([relative_position, variant])
