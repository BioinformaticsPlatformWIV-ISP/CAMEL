import logging
import os
import pandas as pd

from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class GmatsAlgorithm(Tool):

    def __init__(self, camel: Camel) -> None:
        """
        Initializes GmatsAlgorithm.
        :param camel: CAMEL instance
        """
        super().__init__('gMATS Algorithm', '0.1', camel)

    def _execute_tool(self) -> None:
        """
        Executes GmatsAlgorithm.
        :return: text file with gMATS profile
        """

        output_dir = self._parameters['output_directory'].value
        gmats_db = pd.read_csv('/var/lib/.bioit_database/pipelines/neisseria/gMATS_DB.txt', sep='\t')

        for key in self._tool_inputs['TSV']:
            print('Input file : ', key)
            f = open(key.path, "r")

            sample_id = os.path.basename(key.path).replace('.tsv', '').replace('typing-bast-peptide-', '')
            print('sample_ID : ', sample_id)

            output_file = Path(output_dir + '/' + sample_id + '_gMATS.tsv')
            print('Output file : ', output_file)

            # Write output backbone
            with open(output_file, mode='w') as o:
                o.write('Locus' + '\t' +
                        'Allele' + '\t' +
                        '% Identity' + '\t' +
                        'HSP/Locus length' + '\t' +
                        'Allele status' + '\t' +
                        'gMATS status' + '\n')
                f.close()

            # Parse input file
            antigen_info = self.__parse_input_file(key.path)

            # Derive gMATS status for each antigen
            allele_status = []
            for key in antigen_info:
                allele_gmats = self.__derive_allele_gMATS(key, antigen_info[key][1], gmats_db)
                antigen_info[key] += [str(allele_gmats )]
                allele_status.append(allele_gmats)

            # Derive gMATS status of sample
            if 'covered' in allele_status:
                gmats_status = 'covered'
            elif 'not_covered' in allele_status:
                gmats_status = 'not_covered'
            else:
                gmats_status = 'unpredictable'

            # Write output text file
            with open(output_file, mode='a') as o:
                for key in antigen_info:
                    o.write(antigen_info[key][0] + '\t' +
                            antigen_info[key][1] + '\t' +
                            antigen_info[key][2] + '\t' +
                            antigen_info[key][3] + '\t' +
                            antigen_info[key][4] + '\t' +
                            gmats_status + '\n')
                    f.close()

    def __parse_input_file(self, file_path):
        """
        Parses the input file
        :param file_path: Path to input text file
        :return: dictionary with 'Locus', 'Allele', '% identity' for each of the 4 antigens
        """
        allele_info = {}
        with file_path.open() as f:
            for line in f:
                l = line.split('\t')
                if l[0] in ['fHbp_peptide', 'NHBA_peptide', 'NadA_peptide', 'PorA_VR2']:
                    allele_info[l[0]] = l[0:4]
        return allele_info

    def __derive_allele_gMATS(self, locus, allele, gmats_db):
        """
        Derive gMATS prediction for locus/allele
        :param locus: Locus to be assessed
        :param allele: Allele of given locus
        :param gmats_db: Pandas table with the vaccine coverage for the fHbp and NHBA loci
        :return: the gMATS prediction as string
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
    gmats = GmatsAlgorithm(Camel.get_instance())
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







