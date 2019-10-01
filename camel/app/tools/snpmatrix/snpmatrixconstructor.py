import logging
from dataclasses import dataclass
from typing import Dict, List

import os
import vcf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool
# noinspection PyProtectedMember
from vcf.model import _Record as VCFRecord


@dataclass(unsafe_hash=True, frozen=True)
class SNPPosition:
    contig: str
    position: int
    reference_base: chr


class SnpMatrixConstructor(Tool):
    """
    Constructs SNP matrices in FASTA format based on VCF input files.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('SNP Matrix Constructor', '0.1', camel)

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
        print(self._parameters)
        sample_names = [io.value for io in self._tool_inputs['SAMPLE_NAME']]
        nucleotide_by_sample_by_position = self.__get_nucleotides_per_position()
        include_ref = 'include_ref' in self._parameters
        self._tool_outputs['FASTA'] = [self.__generate_matrix(
            sample_names, nucleotide_by_sample_by_position, include_ref)]

    def __get_nucleotides_per_position(self) -> Dict[SNPPosition, Dict[str, str]]:
        """
        Returns a dictionary with the nucleotide for each sample at each variant position.
        :return: Sample_names, nucleotide per position per sample
        """
        nucl_by_position = {}
        for vcf_file, io_sample in zip(self._tool_inputs['VCF'], self._tool_inputs['SAMPLE_NAME']):
            for record in SnpMatrixConstructor.parse_vcf_file(vcf_file.path, 'include_filtered' in self._parameters):
                position = SNPPosition(record.CHROM, record.POS, record.REF)
                if position not in nucl_by_position:
                    nucl_by_position[position] = {}
                nucl_by_position[position][io_sample.value] = str(record.ALT[0])
        logging.info("{} SNP positions found across all samples".format(len(nucl_by_position)))
        return nucl_by_position

    @staticmethod
    def parse_vcf_file(vcf_path: str, include_filtered: bool) -> List[VCFRecord]:
        """
        Parses a single VCF file.
        :param vcf_path: VCF path
        :param include_filtered: If True, filtered variants are included
        :return: List of SNP positions
        """
        vcf_records = []
        with open(vcf_path) as handle:
            for record in list(vcf.Reader(handle))[:50]:
                # Remove non SNP positions
                if not record.is_snp:
                    continue
                # Remove filtered records (unless specified otherwise)
                if (not include_filtered) and len(record.FILTER) > 0:
                    continue
                vcf_records.append(record)
        return vcf_records

    def __generate_matrix(self, sample_names: List[str], nucl_by_pos: Dict[SNPPosition, Dict[str, str]],
                          include_ref: bool = True) -> ToolIOFile:
        """
        Generates a SNP matrix.
        :param sample_names: List of sample names
        :param nucl_by_pos: Nucleotide by sample by SNP position
        :param include_ref: If True, the reference is included in the SNP matrix
        :return: None
        """
        seq_by_sample_name = {name: [] for name in sample_names}
        if include_ref is True:
            seq_by_sample_name['reference'] = []
        for snp_pos, nucl_by_sample in sorted(nucl_by_pos.items(), key=lambda x: x[0].position):
            for name in seq_by_sample_name.keys():
                seq_by_sample_name[name].append(nucl_by_sample.get(name, snp_pos.reference_base))

        output_path = os.path.join(self._folder, 'snp_matrix.fasta')
        seqs = [SeqRecord(Seq(''.join(seq)), name, description='') for name, seq in seq_by_sample_name.items()]
        with open(output_path, 'w') as handle:
            SeqIO.write(seqs, handle, 'fasta')
        return ToolIOFile(output_path)
