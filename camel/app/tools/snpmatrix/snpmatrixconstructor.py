import logging
import os

import vcf

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.snpmatrix.snpposition import SnpPosition
from camel.app.tools.tool import Tool


class SnpMatrixConstructor(Tool):
    """
    Constructs SNP matrices in FASTA format based on VCF input files.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(SnpMatrixConstructor, self).__init__('SNP Matrix Constructor', '0.1', camel)

    def _check_input(self):
        """
        Checks if the tool input is valid.
        :return: None
        """
        if 'VCF' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No VCF input found.")
        super(SnpMatrixConstructor, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        samples_names, snps = self.__parse_vcf_files()
        self._tool_outputs['FASTA'] = [ToolIOFile(self.__generate_matrix(samples_names, snps))]

    def __parse_vcf_files(self):
        """
        Parses the input VCF files.
        :return: Sample_names, snps
        """
        snps = {}
        sample_names = []
        for vcf_file, sample_name in zip(self._tool_inputs['VCF'], self._tool_inputs['SAMPLE_NAME']):
            sample_names.append(sample_name.value)
            vcf_reader = vcf.Reader(open(vcf_file.path))
            for record in vcf_reader:
                if 'INDEL' not in record.INFO:
                    position = SnpPosition(record.CHROM, record.POS, record.REF)
                    if position not in snps:
                        snps[position] = {}
                    snps[position][sample_name.value] = record.ALT[0]
        logging.info("{} SNP positions found across all samples".format(len(snps)))
        return sample_names, snps

    def __generate_matrix(self, samples, snps):
        """
        Generates a SNP matrix.
        :param snps: SNPs
        :return: None
        """
        output_file = os.path.join(self._folder, 'snp_matrix.fasta')
        with open(output_file, 'w') as handle:
            for sample in sorted(samples):
                matrix_row = ''
                for position in sorted(list(snps.keys()), key=lambda p: (p.contig, p.position)):
                    if sample in snps[position]:
                        matrix_row += str(snps[position][sample])
                    else:
                        matrix_row += position.reference_base
                handle.write('>{}\n'.format(sample))
                handle.write(matrix_row)
                handle.write('\n')
        return output_file
