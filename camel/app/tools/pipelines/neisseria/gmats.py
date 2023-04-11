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
        fhbp_db = gmats_db.loc[gmats_db['Locus'] == 'fHbp']
        nhba_db = gmats_db.loc[gmats_db['Locus'] == 'NHBA']

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

            # Prepare lists to collect output
            locus = ['fHbp', 'NHBA', 'NadA', 'PorA_VR2']
            allele = []
            id_pct = []
            coverage = []
            allele_status = []

            # Derive output from input file and gMATS database
            with open(key.path) as f:
                for line in f:
                    l = line.split('\t')
                    print(l)
                    if l[0] == 'fHbp_peptide':
                        allele.append(l[1])
                        id_pct.append(l[2])
                        coverage.append(l[3])
                        if l[1] in ['-', '?']:
                            allele_status.append('unpredictable')
                        elif float(l[1]) in fhbp_db.Allele:
                            allele_status.append(fhbp_db.Status.loc[fhbp_db['Allele'] == float(l[1])].iloc[0])
                        else:
                            allele_status.append('unpredictable')
                    elif l[0] == 'NHBA_peptide':
                        allele.append(l[1])
                        id_pct.append(l[2])
                        coverage.append(l[3])
                        if l[1] in ['-', '?']:
                            allele_status.append('unpredictable')
                        elif float(l[1]) in nhba_db.Allele:
                            print(float(l[1]))
                            print(nhba_db.Status.loc[nhba_db['Allele'] == float(l[1])])
                            allele_status.append(nhba_db.Status.loc[nhba_db['Allele'] == float(l[1])].iloc[0])
                        else:
                            allele_status.append('unpredictable')
                    elif l[0] == 'NadA_peptide':
                        allele.append(l[1])
                        id_pct.append(l[2])
                        coverage.append(l[3])
                        if l[1] in ['-', '?']:
                            allele_status.append('unpredictable')
                        else:
                            allele_status.append('not_covered')
                    elif l[0] == 'PorA_VR2':
                        allele.append(l[1])
                        id_pct.append(l[2])
                        coverage.append(l[3])
                        if l[1] == '4':
                            allele_status.append('covered')
                        elif l[1] in ['-', '?']:
                            allele_status.append('unpredictable')
                        else:
                            allele_status.append('not_covered')

            # Derive gMATS status of sample
            if 'covered' in allele_status:
                gmats_status = ['covered'] * 4
            elif 'not_covered' in allele_status:
                gmats_status = ['not_covered'] * 4
            else:
                gmats_status = ['unpredictable'] * 4

            # Write output text file
            with open(output_file, mode='a') as o:
                for i in range(len(locus)):
                    o.write(locus[i] + '\t' +
                            allele[i] + '\t' +
                            id_pct[i] + '\t' +
                            coverage[i] + '\t' +
                            allele_status[i] + '\t' +
                            gmats_status[i] + '\n')
                    f.close()


#TODO : check peptide integrity via camel
#TODO : Add flag for 'Not covered' alleles with a high identity or coverage length with 'Covered' alleles

if __name__ == '__main__':
    gmats = GmatsAlgorithm(Camel.get_instance())
    logging.info('test')
    gmats.add_input_files({'TSV': [
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-RRS16BD04259.tsv')),
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-S17BD02954.tsv')),
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-S17BD08805.tsv')),
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-S18BD04144.tsv')),
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-S18BD07986.tsv'))]
    })
    logging.info('test2')
    gmats.run()
    logging.info('test3')
    print(gmats.tool_outputs)

#TODO : automatic list of input files







