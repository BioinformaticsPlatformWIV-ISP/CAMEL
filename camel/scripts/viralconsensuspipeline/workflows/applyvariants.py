import dataclasses
import logging
import tempfile
from pathlib import Path
from typing import List, Dict

from Bio import SeqIO

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftoolsconsensus import BcftoolsConsensus
from camel.app.tools.bcftools.bcftoolsindex import BcftoolsIndex
from camel.app.tools.bcftools.bcftoolsview import BcftoolsView


@dataclasses.dataclass
class ApplyVariantsOutput:
    """
    Holder for the output of the 'apply variants' script.
    """
    path_fasta: Path
    informs: List[Dict]


class ApplyVariants(object):
    """
    Wrapper that applies the variants from a VCF file to a FASTA file.
    """

    def __init__(self, dir_: Path) -> None:
        """
        Initializes this workflow.
        :param dir_: Working directory
        :return: None
        """
        self._dir = dir_
        if not self._dir.exists():
            logging.info(f'Creating working directory: {self._dir}')
            self._dir.mkdir(parents=True)
        self._informs = []

    def run(self, fasta_in: Path, vcf_in: Path, name: str, prefix: str = 'consensus', description: str = '') -> \
            ApplyVariantsOutput:
        """
        Runs the variant caller.
        :param fasta_in: Original FASTA file
        :param vcf_in: Input VCF file
        :param name: Dataset name for sequence
        :param prefix: Prefix for the output file
        :param description: Description for the headers
        :return: None
        """
        path_fasta_out = self._dir / f'{prefix}.fasta'
        with tempfile.TemporaryDirectory(prefix='camel_', dir=Camel.get_instance().config['temp_dir']) as dir_temp:
            path_vcf_idx = self.__index_vcf(vcf_in, Path(dir_temp))
            self.__apply_variants(fasta_in, path_vcf_idx, path_fasta_out)
            self.__update_headers(path_fasta_out, name, description)
        return ApplyVariantsOutput(path_fasta_out, self._informs)

    def __index_vcf(self, path_vcf_in: Path, dir_: Path) -> Path:
        """
        Creates an indexed VCF file containing only the sites that passed filtering.
        :param path_vcf_in: Input VCF file
        :param dir_: Working directory
        :return: Path to indexed VCF file
        """
        # Create gzipped VCF file with positions that pass filtering
        path_vcf_out = Path(dir_, f'{path_vcf_in.name}.gz')
        bcftools_view = BcftoolsView(Camel.get_instance())
        bcftools_view.add_input_files({'VCF': [ToolIOFile(path_vcf_in)]})
        bcftools_view.update_parameters(output_type='z', apply_filters='"PASS,."', output_filename=str(path_vcf_out))
        bcftools_view.run(self._dir)
        self._informs.append(bcftools_view)

        # Index VCF file
        bcftools_index = BcftoolsIndex(Camel.get_instance())
        bcftools_index.add_input_files({'VCF_GZ': bcftools_view.tool_outputs['VCF_GZ'] })
        bcftools_index.update_parameters(symlink_input=False)
        bcftools_index.run(self._dir)
        self._informs.append(bcftools_index)

        return path_vcf_out

    def __apply_variants(self, fasta_in: Path, vcf_in: Path, fasta_out: Path) -> None:
        """
        Applies the variants to the target FASTA file.
        :param fasta_in: Input FASTA file
        :param vcf_in: Input VCF file
        :param fasta_out: Output FASTA file
        :return: Path to updated consensus
        """
        bcftools_consensus = BcftoolsConsensus(Camel.get_instance())
        bcftools_consensus.add_input_files({
            'FASTA': [ToolIOFile(Path(fasta_in))],
            'VCF': [ToolIOFile(Path(vcf_in))]
        })
        bcftools_consensus.update_parameters(output_filename=str(fasta_out))
        bcftools_consensus.run(self._dir)
        self._informs.append(bcftools_consensus.informs)

    def __update_headers(self, path_fasta: Path, name: str, description: str = '') -> None:
        """
        Updates the headers of the FASTA file by adding the target description to the headers.
        :param path_fasta: FASTA file to update
        :param name: Dataset name
        :param description: Description to add
        :return: None
        """
        with open(path_fasta) as handle:
            seqs = list(SeqIO.parse(handle, 'fasta'))
        for s in seqs:
            s.id = f"{name}-{s.id.split('-')[-1]}"
            s.description = description
        logging.debug(f'Saving updated header to: {path_fasta}')
        with path_fasta.open('w') as handle:
            SeqIO.write(seqs, handle, 'fasta')
