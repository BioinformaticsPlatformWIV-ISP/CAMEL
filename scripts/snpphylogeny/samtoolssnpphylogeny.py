import logging

import os

from app.components.vcf import vcfutils
from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue
from app.pipeline.pipeline import Pipeline
from app.tools.snpmatrix.snpmatrixconstructor import SnpMatrixConstructor
from resources import YAML_READ_MAPPING_BOWTIE2, YAML_SAMTOOLS_VARIANT_CALLING
from scripts.snpphylogeny.samtools.vcffiltering import VcfFiltering
from scripts.snpphylogeny.snpphylogeny import SnpPhylogeny


class SamtoolsSnpPhylogeny(SnpPhylogeny):
    """
    SNP phylogeny based on the samtools / bcftools SNP calling pipeline.
    """

    def __init__(self):
        """
        Initializes the pipeline script.
        """
        super(SamtoolsSnpPhylogeny, self).__init__('Samtools')
        self._vcf_files_unfiltered = []

    @staticmethod
    def _add_specific_arguments(argument_parser):
        """
        Adds the pipeline specific arguments.
        :param argument_parser: Argument parser
        :return: None
        """
        argument_parser.add_argument('--ploidy', default='haploid', choices=['haploid', 'diploid'])
        argument_parser.add_argument('--calling-method', default='consensus', choices=['consensus', 'multiallelic'])
        argument_parser.add_argument('--min-total-depth', default=10, type=int)
        argument_parser.add_argument('--min-forward-depth', default=1, type=int)
        argument_parser.add_argument('--min-reverse-depth', default=1, type=int)
        argument_parser.add_argument('--min-snp-quality', default=25, type=float)
        argument_parser.add_argument('--min-mapping-quality', default=30, type=int)
        argument_parser.add_argument('--min-distance', default=10, type=int)
        argument_parser.add_argument('--keep-best', action='store_true')
        argument_parser.add_argument('--min-zscore', default=1.96, type=float)
        argument_parser.add_argument('--y-mult', default=10, type=float)

    def _run_snp_calling(self, all_reads):
        """
        Runs the samtools mpileup SNP calling on all input samples.
        :param all_reads: List of PE reads for all samples
        :return: SNP matrix, output files
        """
        for input_reads in all_reads:
            self._bam_files.append(self.__run_read_mapping(input_reads))
        for bam_file in self._bam_files:
            self._vcf_files_unfiltered.append(self.__run_snp_calling(bam_file))
        filtering = VcfFiltering(self._camel, self._destination_path, self._args)
        for sample_name, vcf_file, bam_file in zip(self._sample_names, self._vcf_files_unfiltered, self._bam_files):
            self._vcf_files.append(filtering.apply(sample_name, vcf_file, bam_file))
        self._html.add_filtering_section(self._sample_names, filtering.info)
        snp_matrix = self.__run_snp_matrix_construction()
        return snp_matrix, self.__get_output_files()

    def __run_read_mapping(self, input_reads):
        """
        Runs the read mapping step.
        :param input_reads: Input reads
        :return: BAM file
        """
        pipeline = Pipeline([YAML_READ_MAPPING_BOWTIE2], self._camel, True)
        initial_input = {'FASTA': [ToolIOFile(self._args.reference)]}
        initial_input.update(input_reads)
        pipeline.set_initial_input(initial_input)
        pipeline.run(self._destination_path)
        SamtoolsSnpPhylogeny.__cleanup_read_mapping_files(pipeline)
        return pipeline.outputs['BAM'][0]

    def __run_snp_calling(self, bam_file):
        """
        Runs the SNP calling step.
        :param bam_file: BAM file
        :return: Compressed VCF file 
        """
        pipeline = Pipeline([YAML_SAMTOOLS_VARIANT_CALLING], self._camel, True)
        pipeline.set_initial_input({
            'FASTA': [ToolIOFile(self._args.reference)],
            'BAM': [bam_file]
        })
        job_options = {'skip_variants': 'indels', 'calling_method': self._args.calling_method}
        # Diploid is default setting
        if self._args.ploidy == 'haploid':
            job_options['ploidy'] = '1'
        pipeline.add_job_options({'Variant_calling': job_options})
        pipeline.run(self._destination_path)
        return pipeline.outputs['VCF_GZ'][0]

    def __run_snp_matrix_construction(self):
        """
        Runs the snp matrix construction.
        :return: SNP matrix
        """
        input_files = {
            'SAMPLE_NAME': [ToolIOValue(sample) for sample in self._sample_names],
            'VCF': self._vcf_files
        }
        matrix_constructor = SnpMatrixConstructor(self._camel)
        matrix_constructor.add_input_files(input_files)
        matrix_constructor.run(self._destination_path)
        return matrix_constructor.tool_outputs['FASTA'][0]

    @staticmethod
    def __cleanup_read_mapping_files(pipeline):
        """
        Cleans up the read mapping pipeline intermediate files.
        :param pipeline: Pipeline instance.
        :return: None
        """
        logging.debug("Cleaning up intermediate files from the read mapping pipeline")
        sam_file = pipeline.get_step('Read_mapping').outputs['SAM'][0].path
        logging.debug("Removing {} ({} bytes)".format(sam_file, os.path.getsize(sam_file)))
        os.remove(sam_file)
        bam_unsorted = pipeline.get_step('Sam_to_bam_conversion').outputs['BAM'][0].path
        logging.debug("Removing {} ({} bytes)".format(bam_unsorted, os.path.getsize(bam_unsorted)))
        os.remove(bam_unsorted)

    def __get_output_files(self):
        """
        Returns the output files.
        :return: List with: (Name, Filename, Output files)
        """
        return [
            ('Alignment (BAM)', 'alignment_{}.bam', self._bam_files),
            ('SNPs unfiltered (VCF)', 'snps_{}.vcf.gz', self._vcf_files_unfiltered),
            ('SNPs filtered (VCF)', 'snps_filtered_{}.vcf.gz', self._vcf_files)
        ]

    def _collect_metrics(self):
        """
        Collects the analysis metrics.
        :return: Metrics
        """
        metrics = super(SamtoolsSnpPhylogeny, self)._collect_metrics()
        metrics[0].append('SNPs (Unfiltered)')
        snp_counts = [vcfutils.count_variants(vcf_file.path) for vcf_file in self._vcf_files_unfiltered]
        for row, snps in zip(metrics[1:], snp_counts):
            row.append(snps)
        return metrics

if __name__ == '__main__':
    samtools_snp = SamtoolsSnpPhylogeny()
    samtools_snp.run()
