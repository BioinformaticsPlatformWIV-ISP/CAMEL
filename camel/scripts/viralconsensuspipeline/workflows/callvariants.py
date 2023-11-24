import dataclasses
import logging
from pathlib import Path
from typing import Dict, Any, List

import vcf

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftoolscall import BcftoolsCall
from camel.app.tools.bcftools.bcftoolsmpileup import BcftoolsMpileup
from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
from camel.app.tools.clair3.clair3 import Clair3
from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex


@dataclasses.dataclass
class CallVariantsOutput:
    """
    Holder for the output of the variant calling output.
    """
    path_vcf: Path
    stats: Dict[str, Any]
    informs: List[Dict]


class CallVariants(object):
    """
    Wrapper around the variant callers for the viral consensus pipeline.
    Supported callers: bcftools, Clair3.
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

    def run(self, bam_in: Path, fasta_ref: Path, read_type: str, caller: str, params: Dict[str, Any],
            threads: int = 4) -> CallVariantsOutput:
        """
        Runs the variant caller.
        :param bam_in: Input BAM file
        :param fasta_ref: Reference FASTA file
        :param read_type: Read type
        :param caller: Variant caller (options: GATK, samtools, clair3)
        :param params: Variant calling parameters
        :param threads: Number of threads
        :return: None
        """
        path_ref = self._prepare_ref_genome(fasta_ref)
        if caller == 'bcftools':
            path_vcf = self.call_variants_bcftools(bam_in, path_ref, read_type)
        elif caller == 'clair3':
            path_vcf = self.call_variants_clair3(bam_in, path_ref, read_type, params['model'], threads)
        else:
            raise ValueError(f"Invalid caller: {caller}")
        stats = self._extract_stats(path_vcf)
        return CallVariantsOutput(path_vcf, stats, self._informs)

    def call_variants_bcftools(self, bam_in: Path, path_fasta: Path, read_type: str, max_depth: int = 8_000,
                               threads: int = 4) -> Path:
        """
        Calls variants using bcftools.
        :param path_fasta: Input FASTA file
        :param bam_in: Input BAM file
        :param read_type: Read type
        :param max_depth: Maximum depth
        :param threads: Number of threads
        :return: Path to output VCF file
        """
        # Create pileup
        path_pileup_out = self._dir / 'pileup.vcf.gz'
        bcftools_mpileup = BcftoolsMpileup(Camel.get_instance())
        bcftools_mpileup.update_parameters(
            max_depth=max_depth,
            output_type='z',
            config='illumina' if read_type == 'illumina' else 'ont',
            output_filename=str(path_pileup_out))
        bcftools_mpileup.add_input_files({
            'FASTA': [ToolIOFile(path_fasta)],
            'BAM': [ToolIOFile(bam_in)]
        })
        bcftools_mpileup.run()

        # Call variants
        path_vcf_out = self._dir / 'variants.vcf'
        bcftools_call = BcftoolsCall(Camel.get_instance())
        bcftools_call.add_input_files({'VCF_GZ': [ToolIOFile(path_pileup_out)]})
        bcftools_call.update_parameters(
            output_type='v', output_filename=str(path_vcf_out), ploidy='1', calling_method='consensus',
            variants_only=True, threads=threads)
        bcftools_call.run(self._dir)
        self._informs.extend([bcftools_mpileup.informs, bcftools_call.informs])

        # Save output file
        return path_vcf_out

    def call_variants_clair3(
            self, bam_in: Path, path_fasta: Path, read_type: str, path_model: Path, threads: int) -> Path:
        """
        Calls variants using Clair3.
        :param bam_in: Input BAM file
        :param path_fasta: Input FASTA file
        :param read_type: Read type
        :param path_model: Path to the Clair3 model
        :param threads: Number of threads
        :return: Path to output VCF file
        """
        clair3 = Clair3(Camel.get_instance())
        clair3.add_input_files({'FASTA': [ToolIOFile(path_fasta)], 'BAM': [ToolIOFile(bam_in)]})

        # Check if the model exists
        if not path_model.exists():
            raise FileNotFoundError(f'Clair3 model not found: {path_model}')

        # Options that are always enabled
        platform = 'ilmn' if read_type == 'illumina' else 'ont'
        clair3.update_parameters(
            platform=platform, output_path=str(self._dir), threads=threads, haploid_precise=True, include_ctgs=True,
            no_phasing=True, model_path=str(path_model), chunk_size=5_000)

        # Run Clair3
        clair3.run(self._dir)

        # Convert output VCF.GZ file to VCF format
        path_out = clair3.tool_outputs['VCF'][0].path
        path_out_uncompressed = path_out.parent / path_out.name.replace('.gz', '')
        bcftools_view = BcftoolsView(Camel.get_instance())
        bcftools_view.add_input_files({'VCF_GZ': [ToolIOFile(path_out)]})
        bcftools_view.update_parameters(output_type='v', output_filename=str(path_out_uncompressed))
        bcftools_view.run(self._dir)
        self._informs.append(bcftools_view.informs)
        return path_out_uncompressed

    def _extract_stats(self, path_vcf: Path) -> Dict[str, Any]:
        """
        Extracts variant calling stats by parsing the output VCF file.
        :param path_vcf: Input VCF file
        :return: Variant calling statistics
        """
        variants = list(vcf.Reader(filename=str(path_vcf)))
        return {
            'nb_variants': len(variants),
            'nb_snps': sum(v.is_snp for v in variants),
            'nb_indels': sum(v.is_indel for v in variants)
        }

    def _prepare_ref_genome(self, fasta_in: Path) -> Path:
        """
        Prepares the reference genome by creating a symbolic link in the current directory and creating a samtools
        index.
        :param fasta_in: Input FASTA file
        :return: Path to indexed FASTA file
        """
        # Create symlink
        path_symlink = self._dir / 'ref' / fasta_in.name
        path_symlink.parent.mkdir(exist_ok=True, parents=True)
        if not path_symlink.exists():
            path_symlink.symlink_to(fasta_in)

        # Create FAI index
        samtools_faidx = SamtoolsFastaIndex(Camel.get_instance())
        samtools_faidx.add_input_files({'FASTA': [ToolIOFile(path_symlink)]})
        samtools_faidx.update_parameters(symlink_input=False)
        samtools_faidx.run(self._dir)
        logging.info(f'Reference FASTA file indexed: {path_symlink.name}')
        return path_symlink
