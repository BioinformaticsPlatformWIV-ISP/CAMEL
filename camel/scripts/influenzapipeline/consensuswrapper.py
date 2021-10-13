import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, AbstractSet

import yaml

from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.influenzapipeline import CONFIG_DATA
from camel.scripts.influenzapipeline.snakefile.alignment import OUTPUT_ALIGNMENT_BAM, SNAKEFILE_MAPPING
from camel.scripts.influenzapipeline.snakefile.sequence_extraction import OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE, \
    SNAKEFILE_SEQ_EXTRACTION, OUTPUT_SEQ_EXTRACTION_SELECTVARIANTS, \
    OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE_INDEX_PREFIX


@dataclass
class ConsensusSequenceOutput:
    """
    This dataclass holds the output of the consensus sequence creation workflow.
    """
    fasta_ref: Path
    vcf: Path
    vcf_diffs: Tuple[AbstractSet[str], AbstractSet[str]]
    bam: Path


class ConsensusSequenceWrapper(object):
    """
    This class is used as a wrapper class around the consensus sequence creation Snakemake workflow.
    """

    def __init__(self, working_dir: Path) -> None:
        """
        Initializes the variant calling wrapper.
        :param working_dir: Working directory
        """
        self._working_dir = working_dir
        self._output = None

    @property
    def output(self) -> ConsensusSequenceOutput:
        """
        Returns the filtered VCF file
        :return: VCF file path
        """
        return self._output

    @property
    def input_ref(self) -> Path:
        """
        Returns the path of the input reference that was used
        :return: Path to input reference
        """
        return self._working_dir / 'fasta_ref.io'

    def __create_input(self, fasta_ref_io: Path, fastq_pe: List[ToolIOFile]) -> None:
        """
        Creates the input files for the workflow.
        :param fasta_ref_io: Fasta referencn to use
        :param fastq_pe: PE input files
        :return: None
        """
        fasta_ref = SnakemakeUtils.load_object(fasta_ref_io)[0].path
        SnakemakeUtils.dump_object([ToolIOFile(fasta_ref)], self._working_dir / 'fasta_ref.io')
        SnakemakeUtils.dump_object([ToolIOValue(fasta_ref)], self._working_dir / 'index_genome_prefix.io')
        SnakemakeUtils.dump_object(fastq_pe, self._working_dir / 'fastq_pe.io')

    def __compare_vcf(self, previous_result: Path) -> Tuple[AbstractSet, AbstractSet]:
        """
        Compares the current variants found in the newly created vcf file with those
        from the previous run
        :param previous_result: VCF file from the previous run
        :return: None
        """
        vcf_previous = SnakemakeUtils.load_object(previous_result)[0]
        prev_set = self.__get_variants(vcf_previous)
        vcf_current = SnakemakeUtils.load_object(self._working_dir / OUTPUT_SEQ_EXTRACTION_SELECTVARIANTS)[0]
        curr_set = self.__get_variants(vcf_current)
        return curr_set - prev_set, prev_set - curr_set

    @staticmethod
    def __get_variants(vcf_file: ToolIOFile) -> AbstractSet[Tuple[str, str, str]]:
        """
        Returns the variants found in the VCF file
        :param vcf_file: VCF file to get variants from
        :return: Set of variant tuples (pos, ref, alt)
        """
        with open(vcf_file.path) as inhandle:
            variants = set()
            for line in inhandle:
                if not line.startswith('#'):
                    items = line.split('\t')
                    variant = (items[1], items[3], items[4])  # Position, ref, alt
                    variants.add(variant)
        return variants

    def _get_main_pipeline_config(self):
        """
        Loads the config file from the main pipeline
        :return: Main config dictionary
        """
        with open(CONFIG_DATA) as inhandle:
            return yaml.safe_load(inhandle)

    def _create_config_dictionary(self, aligner: str) -> str:
        """
        Creates the config file for the pipeline
        :param aligner: Aligner to use
        :return: Path to config file
        """
        config_values = {
            'working_dir': str(self._working_dir),
            'aligner': aligner,
            'index_genome_prefix': str(self._working_dir / 'index_genome_prefix.io'),
            'FASTQ_PE': str(self._working_dir / 'fastq_pe.io'),
            'bam': str(self._working_dir / OUTPUT_ALIGNMENT_BAM),
            'fasta_ref': str(self._working_dir / 'fasta_ref.io')
        }

        main_config = self._get_main_pipeline_config()
        if 'iterative_mapping' in main_config['rule_parameters']:
            if aligner == 'bowtie2':
                config_values['rule_parameters'] = {'mapping': main_config['rule_parameters']['iterative_mapping']}

        return SnakePipelineUtils.generate_config_file(config_values, self._working_dir)

    def run_workflow(self, fasta_ref: Path, fastq_pe: List[ToolIOFile], previous_vcf: Path, aligner: str, cores: int = 8) -> None:
        """
        Runs the variant calling workflow.
        :param fasta_ref: Input reference file
        :param fastq_pe: Pre-processed fastq input files
        :param previous_vcf: VCF from the previous consensus sequence creation
        :param aligner: Aligner to use
        :param cores: Number of cores to use for Snakemake
        :return: None
        """
        if not self._working_dir.is_dir():
            self._working_dir.mkdir(parents=True)
        self.__create_input(fasta_ref, fastq_pe)

        # Create config
        config_path = self._create_config_dictionary(aligner)

        # Execute Snakemake for alignment
        output_files = {'BAM': self._working_dir / OUTPUT_ALIGNMENT_BAM}
        SnakePipelineUtils.run_snakemake(SNAKEFILE_MAPPING, config_path, list(output_files.values()), self._working_dir, cores)

        # Execute Snakemake for sequence extraction
        output_files.update({'FASTA': self._working_dir / OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE})
        index_output = self._working_dir / OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE_INDEX_PREFIX
        SnakePipelineUtils.run_snakemake(SNAKEFILE_SEQ_EXTRACTION, config_path, [index_output], self._working_dir, cores)

        self.__set_output(previous_vcf)
        logging.debug(self._output)

    def __set_output(self, previous_vcf: Path) -> None:
        """
        Collects the output of the workflow.
        :param previous_vcf: VCF from the previous consensus sequence creation
        :return: None
        """
        self._output = ConsensusSequenceOutput(
            fasta_ref=self._working_dir / OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE,
            vcf=self._working_dir / OUTPUT_SEQ_EXTRACTION_SELECTVARIANTS,
            vcf_diffs=self.__compare_vcf(previous_vcf),
            bam=self._working_dir / OUTPUT_ALIGNMENT_BAM,
        )
