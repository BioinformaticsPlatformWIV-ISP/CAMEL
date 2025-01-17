import json
from pathlib import Path
from typing import Dict, List, Union

import pandas as pd
import vcf

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.tool import Tool
# noinspection PyProtectedMember
from vcf.model import _Record as VcfRecord

class AMRScreen(Tool):
    """
    Screens the mutations from a VCF file against a database with mutations.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('AMR mutation screen', '0.1', camel)
        self._variant_by_key = None
        self._data_regions = None
        self._ab_short_by_name = None

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'VCF' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Mutation input is required (VCF)")
        if 'VCF_filt' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Filtered mutation input is required (VCF_filt)")
        if 'VCF_lofreq' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Lofreq mutation input is required (VCF_lofreq)")
        if 'DB' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Database input is required (DB)")
        if 'BED' not in self._tool_inputs:
            raise InvalidInputSpecificationError("AMR region input is required (BED)")
        super()._check_input()

    @staticmethod
    def __get_matching_region(position: int, regions: pd.DataFrame) -> Union[Dict, None]:
        """
        Returns the region matching that covers the input position (if available).
        :param position: Genome position
        :param regions: Extracted genomic regions
        :return: The matching region (if available), None otherwise
        """
        for row in regions.to_dict('records'):
            if row['start'] <= position <= row['end']:
                return row
        raise ValueError(f"Position '{position}' does not fall in AMR regions")

    @staticmethod
    def __extracts_mutation_name(row: Dict, full: bool = False) -> str:
        """
        Extracts the mutation name.
        :param row: Input row
        :param full: Full name (including region name)
        :return: Mutation name
        """
        if len(row['associations']) == 0:
            return f'unknown'
        unique_muts = set([x['mutation'] for x in row['associations']])
        if not full:
            return ';'.join(unique_muts) + ('*' if row['passes_filt'] is False else '')
        region_name = row['region']['locus']
        return ';'.join([f'{region_name}_{m}' + ('*' if row['passes_filt'] is False else '') for m in unique_muts])

    def __parse_db_files(self, path_to_db: Path) -> None:
        """
        Parses the TSV files of the database.
        :param path_to_db: Path to database
        :return: None
        """
        # Parse the database with locations
        data_mut_locations = pd.read_table(path_to_db / 'mutation_locations.tsv')
        self._variant_by_key = {
            (r['position'], r['reference_nucleotide'], r['alternative_nucleotide']): r['variant']
            for r in data_mut_locations.to_dict('records')}
        logger.info(f'{len(data_mut_locations):,} mutations parsed')

        # Parse the database with AMR associations
        data_amr_association = pd.read_table(path_to_db / 'amr_associations_all.tsv')
        self._association_by_variant = {}
        for r in data_amr_association.to_dict('records'):
            if r['variant'] not in self._association_by_variant:
                self._association_by_variant[r['variant']] = []
            self._association_by_variant[r['variant']].append(r)
        logger.info(f'{len(data_amr_association):,} AMR associations parsed')

        # Parse antibiotics
        self._data_antibiotics = pd.read_table(path_to_db / 'antibiotics.tsv')
        self._ab_short_by_name = {r['AB']: r['Abbreviation'] for r in self._data_antibiotics.to_dict('records')}

        # Parse AMR regions
        self._data_regions = pd.read_table(self._tool_inputs['BED'][0].path, names=[
            'chr', 'start', 'end', 'locus', 'type', 'abs'])
        self._data_regions['abs_short'] = self._data_regions['abs'].apply(lambda x: ', '.join([
            self._ab_short_by_name[ab] for ab in x.split(', ')]))
        logger.info(f'{len(self._data_regions)} AMR regions parsed')

        # Parse DB version
        try:
            with open(path_to_db / 'VERSION') as handle:
                self._informs['version'] = handle.readline().strip()
        except FileNotFoundError:
            self._informs['version'] = 'n/a'

    @staticmethod
    def __parse_effect(vcf_record: VcfRecord) -> Union[str, None]:
        """
        Parses the mutation effect from the CSQ annotation.
        Note: only extracts it for protein coding regions
        :param vcf_record: Input record
        :return: Mutation effect
        """
        # Check if BCSQ annotation is present
        if 'BCSQ' not in vcf_record.INFO:
            logger.warning(f'BCSQ info missing for: {vcf_record.CHROM}:{vcf_record.POS}')
            return

        # Parse annotation
        parts = vcf_record.INFO['BCSQ'][0].split('|')
        if parts[0].startswith('&'):
            return
        return parts[0]

    def __cross_check_muts_to_db(self, vcf_input: Path, is_lofreq: bool = False) -> List[Dict]:
        """
        Cross-checks the detected mutations to the database.
        :param vcf_input: Input VCF
        :param is_lofreq: is the vcf input from lofreq
        :return: List of detected AMR associations
        """
        mutations_out = []
        with vcf_input.open() as handle:
            variants = list(vcf.Reader(handle))
            logger.info(f"{len(variants):,} variants parsed from '{self._tool_inputs['VCF'][0].path.name}'")

            for vcf_record in variants:
                effect = AMRScreen.__parse_effect(vcf_record)
                region = AMRScreen.__get_matching_region(vcf_record.POS, self._data_regions)
                record = {
                    'alt': ';'.join(str(x) for x in vcf_record.ALT),
                    'effect': effect,
                    'associations': [],
                    'position': vcf_record.POS,
                    'ref': vcf_record.REF,
                    'region': region,
                    'variant_type': vcf_record.var_type,
                    'lofreq': is_lofreq,
                    'passes_filt': False,
                }
                for alt in vcf_record.ALT:
                    # Query the database
                    key = (vcf_record.POS, str(vcf_record.REF), str(alt))
                    if (key not in self._variant_by_key) and effect == 'frameshift':
                        key = f'{region}_LoF'
                    variant = self._variant_by_key.get(key)
                    if variant is None:
                        continue

                    # Add associations
                    for association in self._association_by_variant[variant]:
                        record['associations'].append({
                            'antibiotic': association['drug'],
                            'antibiotic_short': self._ab_short_by_name[association['drug']],
                            'comment': association['comment'] if not pd.isna(association['comment']) else None,
                            'confidence': association['confidence'],
                            'effect': association['effect'],
                            'locus': association['gene'],
                            'mutation': association['mutation'],
                        })
                mutations_out.append(record)
        return mutations_out

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse DB
        self.__parse_db_files(self._tool_inputs['DB'][0].path)

        # Cross-check variants in VCF with DB
        mutations_out = self.__cross_check_muts_to_db(self._tool_inputs['VCF'][0].path)
        logger.info(f"{sum(len(m['associations']) for m in mutations_out):,} AMR associations found")

        # Cross-check variants in VCF with DB separately for Lofreq variants
        mutations_lofreq = self.__cross_check_muts_to_db(self._tool_inputs['VCF_lofreq'][0].path, is_lofreq=True)
        logger.info(f"{sum(len(m['associations']) for m in mutations_lofreq):,} AMR associations found")

        # Check if mutations passed filtering and if they were synonymous
        with open(self._tool_inputs['VCF_filt'][0].path) as handle:
            positions_passing_filt = [record.POS for record in vcf.Reader(handle)]
        mutations_out = [{
            **mut,
            'passes_filt': mut['position'] in positions_passing_filt,
        } for mut in mutations_out]
        logger.info(f"{sum(m['passes_filt'] for m in mutations_out):,}/{len(mutations_out):,} passed filtering")

        # Create text file with mutation positions
        path_bed_out = self.folder / 'mutation_positions.tsv'
        with path_bed_out.open('w') as handle:
            for pos in sorted(set(m['position'] for m in mutations_out)):
                handle.write('\t'.join(['Chromosome', str(pos)]))
                handle.write('\n')
        self._tool_outputs['TSV'] = [ToolIOFile(path_bed_out)]

        # Include the lofreq mutations
        mutations_out.extend(mutations_lofreq)

        # Extract mutation name
        mutations_out = [{
            **row,
            'name': AMRScreen.__extracts_mutation_name(row),
            'name_full': AMRScreen.__extracts_mutation_name(row, full=True),
        } for row in mutations_out]

        # Create JSON output file
        path_json = self.folder / 'amr_screen.json'
        with path_json.open('w') as handle:
            json.dump(mutations_out, handle, indent=2)
        logger.info(f'AMR results stored in: {path_json}')
        self._tool_outputs['JSON'] = [ToolIOFile(path_json)]
