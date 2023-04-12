import logging
import pandas as pd

from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError


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
        if 'output_directory' not in self._parameters:
            raise InvalidInputSpecificationError("Output directory is required.")
        if 'gmats_db' not in self._parameters:
            raise InvalidInputSpecificationError("gMATS database is required")
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("List of input files (TSV) is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the tool.
        :return: 'None'
        """

        output_dir = Path(self._parameters["output_directory"].value)
        gmats_db = pd.read_csv(self._parameters['gmats_db'].value, sep='\t')

        for path_io in self._tool_inputs['TSV']:
            logging.info(f'Processing input file: {path_io}')

            sample_id = path_io.path.stem.replace('typing-bast-peptide-', '')
            logging.info(f'Sample ID: {sample_id}')

            output_file = output_dir / f'{sample_id}_gMATS.tsv'
            logging.info(f'Output file: {output_file}')

            # Write output backbone
            with open(output_file, mode='w') as handle:
                handle.write('\t'.join(['Locus', 'Allele', '% Identity', 'HSP/Locus length', 'Allele status', 'gMATS status', '\n']))

            # Parse input file
            antigen_info = self.__parse_input_file(path_io.path)

            # Derive gMATS status for each antigen
            allele_status = []
            for locus in antigen_info:
                allele_gmats = self.__derive_allele_gmats(locus, antigen_info[locus][1], gmats_db)
                antigen_info[locus] += [str(allele_gmats)]
                allele_status.append(allele_gmats)

            # Derive gMATS status of sample
            if 'covered' in allele_status:
                gmats_status = 'covered'
            elif 'not_covered' in allele_status:
                gmats_status = 'not_covered'
            else:
                gmats_status = 'unpredictable'

            # Write output text file
            with open(output_file, mode='a') as handle:
                for locus in antigen_info:
                    handle.write('\t'.join([antigen_info[locus][0],
                                            str(antigen_info[locus][1]),
                                            str(antigen_info[locus][2]),
                                            str(antigen_info[locus][3]),
                                            antigen_info[locus][4],
                                            gmats_status,
                                            '\n']))
                handle.close()

    def __parse_input_file(self, file_path: Path) -> dict[str, list[str]]:
        """
        Parses the input file.
        :param file_path: Path to input text file
        :return: dictionary with 'Locus', 'Allele', '% identity' for each of the 4 antigens
        """
        input_file = pd.read_table(file_path)

        allele_info = {'fHbp_peptide': input_file.loc[input_file['Locus'] == 'fHbp_peptide'].iloc[0,
                                       0:4].values.flatten().tolist(),
                       'NHBA_peptide': input_file.loc[input_file['Locus'] == 'NHBA_peptide'].iloc[0,
                                       0:4].values.flatten().tolist(),
                       'NadA_peptide': input_file.loc[input_file['Locus'] == 'NadA_peptide'].iloc[0,
                                       0:4].values.flatten().tolist(),
                       'PorA_VR2': input_file.loc[input_file['Locus'] == 'PorA_VR2'].iloc[0,
                                   0:4].values.flatten().tolist()}

        return allele_info

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
            allele_gmats = 'unpredictable'
        elif locus == 'fHbp_peptide':
            if float(allele) in fhbp_db.Allele:
                allele_gmats = fhbp_db.Status.loc[fhbp_db['Allele'] == float(allele)].iloc[0]
            else:
                allele_gmats = 'unpredictable'
        elif locus == 'NHBA_peptide':
            if float(allele) in nhba_db.Allele:
                allele_gmats = nhba_db.Status.loc[nhba_db['Allele'] == float(allele)].iloc[0]
            else:
                allele_gmats = 'unpredictable'
        elif locus == 'NadA_peptide':
            allele_gmats = 'not_covered'
        elif locus == 'PorA_VR2':
            allele_gmats = 'covered' if allele == '4' else 'not_covered'
        return allele_gmats

#TODO : check peptide integrity via camel
#TODO : Add flag for 'Not covered' alleles with a high identity or coverage length with 'Covered' alleles


if __name__ == '__main__':
    gmats = GMats(Camel.get_instance())
    gmats.add_input_files({'TSV': [
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-RRS16BD04259.tsv')),
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-S17BD02954.tsv')),
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-S17BD08805.tsv')),
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-S18BD04144.tsv')),
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-S18BD07986.tsv'))]
    })
    gmats.run()
    logging.info('Done')

#TODO : automatically list input files