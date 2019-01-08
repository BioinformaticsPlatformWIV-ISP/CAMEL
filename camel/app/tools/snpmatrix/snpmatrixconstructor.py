import logging
from typing import Dict, List, Tuple

import os
import vcf

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.snpmatrix.snpposition import SnpPosition
from camel.app.tools.tool import Tool


class SnpMatrixConstructor(Tool):
    """
    Constructs SNP matrices in FASTA format based on VCF input files.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(SnpMatrixConstructor, self).__init__('SNP Matrix Constructor', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the tool input is valid.
        :return: None
        """
        if 'VCF' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No VCF input found.")
        super(SnpMatrixConstructor, self)._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        samples_names, snps = self.__parse_vcf_files()
        include_ref = self._parameters['include_ref'].as_boolean()
        self._tool_outputs['FASTA'] = [self.__generate_matrix(samples_names, snps, include_ref)]

    def __parse_vcf_files(self) -> Tuple[List[str], Dict[SnpPosition, Dict[str, str]]]:
        """
        Parses the input VCF files.
        :return: Sample_names, snps
        """
        nucl_by_position = {}
        sample_names = []
        for vcf_file, sample_name in zip(self._tool_inputs['VCF'], self._tool_inputs['SAMPLE_NAME']):
            sample_names.append(sample_name.value)
            with open(vcf_file.path) as handle:
                vcf_reader = vcf.Reader(handle)
                for record in vcf_reader:
                    if 'INDEL' in record.INFO:
                        continue
                    position = SnpPosition(record.CHROM, record.POS, record.REF)
                    if position not in nucl_by_position:
                        nucl_by_position[position] = {}
                    nucl_by_position[position][sample_name.value] = str(record.ALT[0])
        logging.info("{} SNP positions found across all samples".format(len(nucl_by_position)))
        return sample_names, nucl_by_position

    def __generate_matrix(self, sample_names: List[str], nucl_by_pos: Dict[SnpPosition, Dict[str, str]],
                          include_ref: bool = True) -> ToolIOFile:
        """
        Generates a SNP matrix.
        :param sample_names: List of sample names
        :param nucl_by_pos: Nucleotide per sample by SNP position
        :param include_ref: If True, the reference is included in the SNP matrix
        :return: None
        """
        seq_by_sample_name = {name: [] for name in sample_names}
        if include_ref is True:
            seq_by_sample_name['reference'] = []
        for snp_pos, nucl_by_sample in sorted(nucl_by_pos.items()):
            for name in seq_by_sample_name.keys():
                seq_by_sample_name[name].append(nucl_by_sample.get(name, snp_pos.reference_base))

        output_path = os.path.join(self._folder, 'snp_matrix.fasta')
        with open(output_path, 'w') as handle:
            for name, sequence in seq_by_sample_name.items():
                handle.write(f'>{name}')
                handle.write('\n')
                handle.write(''.join(sequence))
                handle.write('\n')
        return ToolIOFile(output_path)
