import logging
from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class GMats(Tool):
    """
    gMATS is the genetic Meningococcal Antigen Typing System (gMATS) for Neisseria meningitidis
    used to predict the Bexsero vaccine effectiveness.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the GMats tool.
        :param camel: CAMEL instance
        """
        super().__init__('gMATS Algorithm', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'DB' not in self._tool_inputs:
            raise InvalidInputSpecificationError("gMATS database is required")
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("List of input files (TSV) is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the tool.
        :return: 'None'
        """
        gmats_db = pd.read_table(self._tool_inputs['DB'][0].path)

        # antigen_info = self._parse_input_file()
        # self._determine_gmats_alelle(antigen_info, db)
        # self._determine_gmats_status(antigen_info, db)
        # self._create_output_file(antigen_info)

        # Parse input file
        tsv_in = self._tool_inputs['TSV'][0].path
        sample_id = tsv_in.stem.replace('typing-bast-peptide-', '')
        logging.info(f'Sample ID: {sample_id}')
        antigen_info = self.__parse_input_file(tsv_in)

        # Determine gMATS allele
        antigen_info['gmats_allele'] = antigen_info.apply(
            lambda x: self.__derive_allele_gmats(x['Locus'], x['Allele'], gmats_db), axis=1)

        # Derive gMATS status of sample
        if 'covered' in list(antigen_info['gmats_allele']):
            antigen_info['gMATS status'] = 'covered'
        elif 'not_covered' in list(antigen_info['gmats_allele']):
            antigen_info['gMATS status'] = 'not_covered'
        else:
            antigen_info['gMATS status'] = 'unpredictable'

        # Get the output directory
        if 'output_directory' in self._parameters:
            output_dir = Path(self._parameters["output_directory"].value)
        else:
            output_dir = self.folder

        # Create output file
        output_file = output_dir / f'{sample_id}_gMATS.tsv'
        antigen_info.drop(['Type'], axis=1, inplace=True)
        antigen_info.to_csv(output_file, sep='\t', index=False)
        logging.info(f'Output file created: {output_file}')
        self._tool_outputs['TSV'] = [ToolIOFile(output_file)]

    def __parse_input_file(self, file_path: Path) -> pd.DataFrame:
        """
        Parses the input file.
        :param file_path: Path to input text file
        :return: DataFrame with allele info
        """
        input_file = pd.read_table(file_path, dtype=str)
        input_file = input_file[input_file['Locus'] != 'PorA_VR1']
        return input_file

    def __derive_allele_gmats(self, locus: str, allele: str, gmats_db: pd.DataFrame) -> str:
        """
        Derive gMATS prediction for locus/allele.
        :param locus: Locus to be assessed
        :param allele: Specific allele to be assessed
        :param gmats_db: Pandas table containing the vaccine coverage for the fHbp and NHBA loci
        :return: the gMATS prediction
        """
        fhbp_db = gmats_db.loc[gmats_db['Locus'] == 'fHbp']
        nhba_db = gmats_db.loc[gmats_db['Locus'] == 'NHBA']

        if allele in ['-', '?']:
            return 'unpredictable'

        # fHbp
        if locus == 'fHbp_peptide':
            if float(allele) in fhbp_db.Allele:
                return fhbp_db.Status.loc[fhbp_db['Allele'] == float(allele)].iloc[0]
            else:
                return 'unpredictable'

        # NHBA
        if locus == 'NHBA_peptide':
            if float(allele) in nhba_db.Allele:
                return nhba_db.Status.loc[nhba_db['Allele'] == float(allele)].iloc[0]
            else:
                return 'unpredictable'

        # PorA_VR2
        if locus == 'PorA_VR2':
            return 'covered' if allele == '4' else 'not_covered'

        # Other loci
        if locus in ('NadA_peptide', 'PorA_VR1'):
            return 'not_covered'
        raise ValueError(f'Cannot determine gmats allele (locus={locus}, allele={allele})')

        # TODO: check peptide integrity via camel
        # TODO: Add flag for 'Not covered' alleles with a high identity or coverage length with 'Covered' alleles
