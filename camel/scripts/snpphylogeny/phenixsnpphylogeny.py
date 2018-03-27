import os

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.phenix.phenixpreparereference import PhenixPrepareReference
from camel.app.tools.phenix.phenixsnppipeline import PhenixSnpPipeline
from camel.app.tools.phenix.phenixvcfs2fasta import PhenixVcfs2Fasta
from camel.scripts.snpphylogeny.snpphylogeny import SnpPhylogeny


class PhenixSnpPhylogeny(SnpPhylogeny):
    """
    SNP phylogeny based on the PHEnix SNP calling pipeline.
    """

    THREADS = 8

    def __init__(self):
        """
        Initializes the pipeline script.
        """
        super(PhenixSnpPhylogeny, self).__init__('PHEnix')

    @staticmethod
    def _add_specific_arguments(argument_parser):
        """
        Adds the pipeline specific arguments.
        :param argument_parser: Argument parser
        :return: None
        """
        argument_parser.add_argument('--mapper', choices=['bwa', 'bowtie2'], default='bwa')
        argument_parser.add_argument('--variant-caller', choices=['gatk', 'mpileup'], default='gatk')
        argument_parser.add_argument('--min-snp-quality', default=25, type=int)
        argument_parser.add_argument('--min-mapping-quality', default=30, type=int)
        argument_parser.add_argument('--min-snp-depth', default=10, type=int)
        argument_parser.add_argument('--ploidy')
        argument_parser.add_argument('--column-gaps', type=float)
        argument_parser.add_argument('--column-ns', type=float)

    def _run_snp_calling(self, reads):
        """
        Runs the SNP calling.
        :param reads: List of read inputs for the SNP calling
        :return: SNP matrix, output files
        """
        reference = self.__prepare_reference()
        filter_string = 'min_depth:{},mq_score:{},qual_score:{}'.format(
            self._args.min_snp_depth, self._args.min_mapping_quality, self._args.min_snp_quality)

        for sample_name, sample_reads in zip(self._sample_names, reads):
            snp_pipeline = PhenixSnpPipeline(self._camel)
            snp_pipeline.add_input_files({
                'FASTA': [reference],
                'FASTQ_PE': sample_reads['FASTQ_PE'],
                'SAMPLE_NAME': [ToolIOValue(sample_name)]
            })

            if self._args.mapper == 'bowtie2':
                snp_pipeline.update_parameters(mapper_options='"--threads {}"'.format(PhenixSnpPhylogeny.THREADS))
            else:
                snp_pipeline.update_parameters(mapper_options='"-t {}"'.format(PhenixSnpPhylogeny.THREADS))

            if self._args.variant_caller == 'gatk':
                variant_options = '"-ploidy {} --num_threads {}"'.format(
                    {'haploid': 1, 'diploid': 2}.get(self._args.ploidy, 2), PhenixSnpPhylogeny.THREADS)
            else:
                variant_options = False
            snp_pipeline.update_parameters(mapper=self._args.mapper, variant=self._args.variant_caller,
                                           variant_options=variant_options, filters=filter_string)
            sample_dir = os.path.join(self._destination_path, FileSystemHelper.make_valid(sample_name))
            os.mkdir(sample_dir)
            snp_pipeline.run(sample_dir)

            self._bam_files.append(snp_pipeline.tool_outputs['BAM'][0])
            if 'VCF_Filt' in snp_pipeline.tool_outputs:
                vcf_file = snp_pipeline.tool_outputs['VCF_Filt'][0]
            else:
                vcf_file = snp_pipeline.tool_outputs['VCF'][0]
            self._vcf_files.append(vcf_file)
        output_files = self.__get_output_files()
        snp_matrix = self.__generate_snp_matrix()
        return snp_matrix, output_files

    def __prepare_reference(self):
        """
        Prepares the reference.
        :return: None
        """
        prepare_reference = PhenixPrepareReference(self._camel)
        prepare_reference.add_input_files({'FASTA': [ToolIOFile(self._args.reference)]})
        prepare_reference.update_parameters(mapper=self._args.mapper, variant=self._args.variant_caller,
                                            reference_name=self._args.reference_name)
        prepare_reference.run(self._destination_path)
        return prepare_reference.tool_outputs['FASTA'][0]

    def __generate_snp_matrix(self):
        """
        Generates the SNP matrix.
        :return: None
        """
        vcf_to_fasta = PhenixVcfs2Fasta(self._camel)
        vcf_to_fasta.add_input_files({'VCF': self._vcf_files})
        vcf_to_fasta.update_parameters(column_gaps=self._args.column_gaps, column_Ns=self._args.column_ns)
        vcf_to_fasta.run(self._destination_path)
        return vcf_to_fasta.tool_outputs['FASTA'][0]

    def __get_output_files(self):
        """
        Returns the output files.
        :return: List with (Name, Filename template, Output files)
        """
        return [
            ('Alignment (BAM)', 'alignment_{}.bam', self._bam_files),
            ('SNPs (VCF)', 'snps_{}.vcf.gz', self._vcf_files)
        ]


if __name__ == '__main__':
    phenix_snp = PhenixSnpPhylogeny()
    phenix_snp.run()
