import argparse
import re
import shutil
from pathlib import Path
from typing import Optional, Sequence, Union, List, Dict

import pandas as pd

from camel.app.camel import Camel
from camel.app.loggers import logger


class UpdateAMRDB(object):
    """
    Class to update the Mycobacterium AMR database.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the updater.
        :param args: Command line arguments
        :return: None
        """
        self._args = UpdateAMRDB._parse_arguments(args)
        self._data_mut_locations = None
        self._data_catalogue = None
        self._data_antibiotics = None
        self._gff_by_gene_name = {}
        self._gff_by_locus_tag = {}

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments (optional)
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        # Input
        parser.add_argument(
            '--tsv-coords', required=True, type=Path,
            help='TSV file with the genomic coordinates of the mutations')
        parser.add_argument(
            '--tsv-catalogue', required=True, type=Path,
            help='TSV file with the mutation catalogue')
        parser.add_argument(
            '--tsv-abs', required=True, type=Path,
            help='TSV file with antibiotics')
        parser.add_argument(
            '--ref-gff', required=True, type=Path,
            help='GFF3 file with the ref. genome annotations, used to extract regions')
        parser.add_argument('--version', type=str, help='DB version', required=True)

        # Other parameters
        parser.add_argument('--dir-out', type=Path, help='Output directory')
        return parser.parse_args(args)

    def __parse_input_tsv_files(self) -> None:
        """
        Parses the input TSV files.
        :return: None
        """
        # Parse the database with locations
        self._data_mut_locations = pd.read_table(self._args.tsv_coords)
        logger.info(f'{len(self._data_mut_locations):,} genomic coordinates parsed')

        # Parse the database with the AMR catalogue
        self._data_catalogue = pd.read_table(self._args.tsv_catalogue)
        logger.info(f'{len(self._data_catalogue):,} AMR associations parsed')

        # Parse antibiotics
        self._data_antibiotics = pd.read_table(self._args.tsv_abs)
        logger.info(f'{len(self._data_antibiotics):,} Antibiotics parsed')

    def __parse_input_gff(self) -> None:
        """
        Parses the input GFF reference genome annotation.
        :return: None
        """
        columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        data_gff = pd.read_table(self._args.ref_gff, names=columns, comment='#', header=None, sep='\t')
        logger.info(f'{len(data_gff):,} features parsed from input GFF file')

        # Retain only 'gene' features
        data_gff = data_gff[data_gff['type'] == 'gene']
        logger.info(f"{len(data_gff):,} features marked as 'gene'")

        # Store in dictionaries for retrieval
        data_gff['gene_name'] = data_gff['attributes'].apply(lambda x: re.search('gene=(.*?);', x).group(1))
        data_gff['locus_tag'] = data_gff['attributes'].apply(lambda x: re.search('locus_tag=(.*?);', x).group(1))
        self._gff_by_gene_name = {r['gene_name']: r for r in data_gff.to_dict('records')}
        self._gff_by_locus_tag = {r['locus_tag']: r for r in data_gff.to_dict('records')}

    def __extract_gene_regions(self) -> pd.DataFrame:
        """
        Extracts the regions with the AMR-associated genes from the parsed GFF information.
        :return: Data frame with region data
        """
        records_out = []
        for gene_name, data in self._data_catalogue.groupby('gene'):
            logger.debug(f"Processing '{gene_name}' ({len(data):,} associations)")

            # Get GFF information for gene
            try:
                row_gff = self._gff_by_gene_name[gene_name]
            except KeyError:
                row_gff = self._gff_by_locus_tag[gene_name]

            # Store record
            records_out.append({
                'chr': 'Chromosome',
                'start': row_gff['start'],
                'end': row_gff['end'],
                'gene': gene_name,
                'type': re.search('gene_biotype=(.*?);', row_gff['attributes']).group(1),
                'abs': ', '.join(sorted(data['drug'].unique()))
            })
        return pd.DataFrame(records_out)

    def __extract_prom_regions(self) -> pd.DataFrame:
        """
        Extracts the regions with the promoters of the AMR-associated genes that contain mutations from the catalogue.
        :return: Data frame with region data
        """
        # Merge data frames (note that a single variant can match multiple positions)
        data_merged = pd.merge(self._data_catalogue, self._data_mut_locations, how='left', on='variant')

        # Identify promotor mutations
        data_merged['is_promotor'] = data_merged['mutation'].apply(lambda x: x.startswith('c.-') or x.startswith('n.-'))

        # Add regions
        records_out = []
        for gene, data_mutations in data_merged[data_merged['is_promotor']].groupby('gene'):
            p_start, p_end = (int(data_mutations['position'].min()), int(data_mutations['position'].max()))
            logger.debug(f'{len(data_mutations):,} promotor mutations for {gene}')
            records_out.append({
                'chr': 'Chromosome',
                'start': p_start,
                'end': p_end,
                'gene': f'{gene}_prom',
                'type': 'promotor',
                'abs': ', '.join(sorted(data_mutations['drug'].unique()))
            })
        return pd.DataFrame(records_out)

    @staticmethod
    def __get_matching_region(position: int, regions: List[Dict]) -> Union[str, None]:
        """
        Returns the region matching that covers the input position (if available).
        :param position: Genome position
        :param regions: Extracted genomic regions
        :return: The matching region (if available), None otherwise
        """
        for row in regions:
            if row['start'] <= position <= row['end']:
                return row['gene']
        return None

    def run(self) -> None:
        """
        Runs the updater.
        :return: None
        """
        # Parse input file
        self.__parse_input_tsv_files()
        self.__parse_input_gff()

        # Construct BED file with target gene regions
        data_regions_genes = self.__extract_gene_regions()
        logger.info(f'{len(data_regions_genes):,} gene regions extracted')

        # As the promoter mutations fall outside the gene regions, they need to be extracted separately
        data_regions_prom = self.__extract_prom_regions()
        logger.info(f'{len(data_regions_genes):,} promotor regions extracted')

        # Export merged BED file
        data_regions_all = pd.concat([data_regions_genes, data_regions_prom])
        data_regions_all.sort_values(by='start', inplace=True)
        data_regions_all.to_csv(self._args.dir_out / 'amr_regions.bed', sep='\t', header=False, index=False)
        logger.info('BED file created')

        # Extract mutations that fall outside the extracted regions
        regions_list = data_regions_all.to_dict('records')
        self._data_mut_locations['region'] = self._data_mut_locations['position'].apply(
            lambda x: UpdateAMRDB.__get_matching_region(x, regions_list))
        mutations_not_covered = self._data_mut_locations[pd.isna(self._data_mut_locations['region'])]
        logger.info(f"{len(mutations_not_covered):,}/{len(self._data_mut_locations):,} mutations not covered")
        mutations_not_covered.to_csv(self._args.dir_out / '_not_covered.tsv', sep='\t', index=False)

        # Verify that the antibiotics in the AMR catalogue are present in the TSV file
        short_by_full = {ab['AB']: ab for ab in self._data_antibiotics.to_dict('records')}
        self._data_catalogue['drug'].apply(lambda x: short_by_full[x])
        logger.info(f'All antibiotics from AMR catalogue are present in the input TSV file')

        # Copy the TSV files to the output folder
        shutil.copyfile(self._args.tsv_coords, self._args.dir_out / 'mutation_locations.tsv')
        shutil.copyfile(self._args.tsv_catalogue, self._args.dir_out / 'amr_associations.tsv')
        shutil.copyfile(self._args.tsv_abs, self._args.dir_out / 'antibiotics.tsv')

        # Create a text file with the DB version
        with open(self._args.dir_out / 'VERSION', 'w') as handle:
            handle.write(self._args.version)
            handle.write('\n')
        logger.info(f'Version file created')


if __name__ == '__main__':
    Camel.get_instance()
    updater = UpdateAMRDB()
    updater.run()
