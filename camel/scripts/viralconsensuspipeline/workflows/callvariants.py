import dataclasses
import logging
import re
from pathlib import Path
from typing import Dict, Any, List

import vcf

from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.clair3.clair3 import Clair3


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

    def run(self, bam_in: Path, fasta_ref: Path, input_type: str, caller: str, params: Dict[str, Any],
            threads: int = 4) -> CallVariantsOutput:
        """
        Runs the variant caller.
        :param bam_in: Input BAM file
        :param fasta_ref: Reference FASTA file
        :param input_type: Read type
        :param caller: Variant caller (options: GATK, samtools, clair3)
        :param params: Variant calling parameters
        :param threads: Number of threads
        :return: None
        """
        path_ref = self._prepare_ref_genome(fasta_ref)
        if caller == 'bcftools':
            path_vcf = self.call_variants_bcftools(bam_in, path_ref, input_type)
        elif caller == 'clair3':
            path_vcf = self.call_variants_clair3(bam_in, path_ref, input_type, params['model'], threads)
        else:
            raise ValueError(f"Invalid caller: {caller}")
        stats = self._extract_stats(path_vcf)
        return CallVariantsOutput(path_vcf, stats, self._informs)

    def call_variants_bcftools(self, bam_in: Path, path_fasta: Path, input_type: str, max_depth: int = 8_000) -> Path:
        """
        Calls variants using bcftools.
        :param path_fasta: Input FASTA file
        :param bam_in: Input BAM file
        :param input_type: Read type
        :param max_depth: Maximum depth
        :return: Output VCF file
        """
        # Create working directory
        path_vcf_out = self._dir / 'variants.vcf'

        # Config
        config = 'illumina' if input_type == 'illumina' else 'ont'

        # Collect input file
        command = Command(' '.join([
            'module load bcftools/1.17;',
            'bcftools mpileup', f'-d {max_depth}', '--output-type v', f'--config {config}',
            f'--fasta-ref {path_fasta}', str(bam_in), '|',
            f'bcftools call --output {path_vcf_out} --output-type v --ploidy 1 --consensus-caller --variants-only'
        ]))
        command.run(self._dir)
        if not command.returncode == 0:
            raise RuntimeError(command.stderr)
        self._informs.append({'_name': 'bcftools mpileup 1.17', '_version': '1.17', '_command': command.command})

        # Save output file
        return path_vcf_out

    def call_variants_clair3(
            self, bam_in: Path, path_fasta: Path, input_type: str, path_model: Path, threads: int) -> Path:
        """
        Calls variants using Clair3.
        :param bam_in: Input BAM file
        :param path_fasta: Input FASTA file
        :param input_type: Read type
        :param path_model: Path to the Clair3 model
        :param threads: Number of threads
        """
        clair3 = Clair3(Camel.get_instance())
        clair3.add_input_files({'FASTA': [ToolIOFile(path_fasta)], 'BAM': [ToolIOFile(bam_in)]})

        # Check if the model exists
        if not path_model.exists():
            raise FileNotFoundError(f'Clair3 model not found: {path_model}')

        # Options that are always enabled
        platform = 'ilmn' if input_type == 'illumina' else 'ont'
        clair3.update_parameters(
            platform=platform, output_path=str(self._dir), threads=threads, haploid_precise=True, include_ctgs=True,
            no_phasing=True, model_path=str(path_model), chunk_size=5_000)

        # Run Clair3
        clair3.run(self._dir)

        # Convert output VCF.GZ file to VCF format
        path_out = clair3.tool_outputs['VCF'][0].path
        path_out_uncompressed = path_out.parent / path_out.name.replace('.gz', '')
        command = Command(f'module load bcftools/1.17; bcftools view {path_out} > {path_out_uncompressed};')
        command.run(self._dir)
        self._informs.append(clair3.informs)
        path_out_uncompressed = self.__filter_problematic_positions(path_out_uncompressed)
        return path_out_uncompressed

    def __filter_problematic_positions(self, vcf_in: Path) -> Path:
        """
        Filters the problematic positions reported by Clair3.
        :param vcf_in: Input file in VCF format
        :return: VCF file with problematic positions removed
        """
        format_ok = True
        try:
            with open(vcf_in) as handle:
                list(vcf.Reader(handle))
        except ValueError as err:
            format_ok = False
            logging.warning(f'Could not parse VCF generated by Clair3: {err}')

        # Parsing is OK -> return original file
        if format_ok:
            return vcf_in

        # Parsing is not OK -> remove problematic lines
        lines_out = []
        with open(vcf_in) as handle:
            for line in handle.readlines():
                # Header lines
                if line.startswith('#'):
                    lines_out.append(line)
                    continue
                parts = line.strip().split('\t')
                m = re.search(':\d+\.\d+,\d+$', parts[-1])
                if not m:
                    lines_out.append(line)
                    continue
                logger.warning(f'Skipping variant at position {parts[1]} - incorrect format')

        # Create output file
        path_out = vcf_in.parent / vcf_in.name.replace('.vcf', '-fixed.vcf')
        with open(path_out, 'w') as handle:
            handle.write(''.join(lines_out))
        return path_out

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
        Prepares the reference genome by creating a symbolic link in the  current directory and creating a samtools
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
        command = Command(f'module load samtools/1.17; samtools faidx {path_symlink};')
        command.run(self._dir)
        if not command.returncode == 0:
            raise RuntimeError(f'Error creating FASTA index: {command.stderr}')
        logging.info(f'Reference FASTA file indexed: {path_symlink.name}')
        return path_symlink
