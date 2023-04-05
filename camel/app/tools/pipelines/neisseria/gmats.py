import logging
import os

from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class GmatsAlgorithm(Tool):

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('gMATS Algorithm', '0.1', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: concatenated tsv files
        """

        output_dir = self._parameters['output_directory'].value

        for key in self._tool_inputs['TSV']:
            print(key)
            f = open(key.path, "r")

            sample_ID = os.path.basename(key.path).replace('.tsv', '').replace('typing-bast-peptide-', '')
            print('Sample_ID : ', sample_ID)

            output_file = Path(output_dir + '/' + sample_ID + '_gMATS.tsv')
            print('Output file : ', output_file)

            with open(output_file, mode = 'w') as o:
                o.write('Locus' + '\t' +
                        'Allele' + '\t' +
                        '% Identity' + '\t' +
                        'HSP/Locus length' + '\t' +
                        'Allele status' + '\t' +
                        'gMATS status' + '\n')
                f.close()

            locus = ['fHbp', 'NHBA', 'NadA', 'PorA_VR2']
            allele = []
            id_pct = []
            coverage = []
            allele_status = []
            gmats_status = ['-', '-', '-', '-']

            with open(key.path) as f:
                for line in f:
                    l = line.split('\t')
                    print(l)
                    if l[0] == 'fHbp_peptide':
                        allele.append(l[1])
                        id_pct.append(l[2])
                        coverage.append(l[3])
                        allele_status.append('')
                    elif l[0] == 'NHBA_peptide':
                        allele.append(l[1])
                        id_pct.append(l[2])
                        coverage.append(l[3])
                        allele_status.append('')
                    elif l[0] == 'NadA_peptide':
                        allele.append(l[1])
                        id_pct.append(l[2])
                        coverage.append(l[3])
                        allele_status.append('Not covered')
                    elif l[0] == 'PorA_VR2':
                        allele.append(l[1])
                        id_pct.append(l[2])
                        coverage.append(l[3])
                        if l[1] == '4':
                            allele_status.append('Covered')
                        else:
                            allele_status.append('Not covered')

            print('Loci :' + str(locus))
            print('Alleles :' + str(allele))
            print('ID_pct :' + str(id_pct))
            print('Coverage :' + str(coverage))
            print('Allele status' + str(allele_status))
            print('gMATS status' + str(gmats_status))

            with open(output_file, mode = 'a') as o:
                for i in range(len(locus)):
                    print(i)
                    o.write(locus[i] + '\t' +
                            allele[i] + '\t' +
                            id_pct[i] + '\t' +
                            coverage[i] + '\t' +
                            allele_status[i] + '\t' +
                            gmats_status[i] + '\n')
                    f.close()

#TODO : infer allele status -> best way to proceed?
#TODO : infer sample gMATS status


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
