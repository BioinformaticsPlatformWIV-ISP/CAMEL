import argparse

from app.camel import Camel
from app.io.tooliofile import ToolIOFile
from app.pipeline.pipeline import Pipeline
from resources import YAML_SAMTOOLS_VARIANT_CALLING


class SamtoolsVariantCalling(object):
    """
    Class to run samtools variant calling using CAMEL.
    """

    def __init__(self):
        """
        Initializes the main script.
        """
        self._args = SamtoolsVariantCalling._parse_arguments()
        self._camel = Camel()
        self._pipeline = Pipeline([YAML_SAMTOOLS_VARIANT_CALLING], self._camel, True)

    @staticmethod
    def _parse_arguments():
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--bam', required=True)
        argument_parser.add_argument('--reference', required=True)
        argument_parser.add_argument('--output', required=True)
        argument_parser.add_argument('--ploidy', choices=['GRCh37', 'GRCh38', 'X', 'Y', '1'])
        argument_parser.add_argument('--skip-variants', choices=['snps', 'indels'])
        argument_parser.add_argument('--count-orphans', action='store_true')
        argument_parser.add_argument('--disable-baq', action='store_true')
        argument_parser.add_argument('--minimal-mq', type=int)
        argument_parser.add_argument('--minimal-bq', type=int)
        argument_parser.add_argument('--calling-method', choices=('consensus', 'multiallelic'))
        argument_parser.add_argument('--output-all-sites', action='store_true')
        argument_parser.add_argument('--mutation-rate')
        return argument_parser.parse_args()

    def call_variants(self):
        """
        Call SNPs on the given BAM file.
        :return: None
        """
        self.__set_initial_input()
        self.__update_pileup_parameters()
        self.__update_variant_calling_parameters()
        self._pipeline.run('.')
        self.__save_output_file()

    def __set_initial_input(self):
        """
        Sets the initial input of the pipeline.
        :return: None
        """
        self._pipeline.set_initial_input({
            'FASTA': [ToolIOFile(self._args.reference)],
            'BAM': [ToolIOFile(self._args.bam)]
        })

    def __update_pileup_parameters(self):
        """
        Updates the pileup parameters.
        :return: None
        """
        if self._args.minimal_mq != 0:
            self._pipeline.add_job_options({'Pileup': {'min_mapping_quality': self._args.minimal_mq}})
        if self._args.minimal_bq != 0:
            self._pipeline.add_job_options({'Pileup': {'min_base_quality': self._args.minimal_bq}})
        if self._args.count_orphans is True:
            self._pipeline.add_job_options({'Pileup': {'count_orphans': True}})
        if self._args.disable_baq is True:
            self._pipeline.add_job_options({'Pileup': {'disable_baq': True}})

    def __update_variant_calling_parameters(self):
        """
        Updates the variant calling parameters.
        :return: None
        """
        if self._args.ploidy:
            self._pipeline.add_job_options({'Variant_calling': {'ploidy': self._args.ploidy}})
        if self._args.skip_variants:
            self._pipeline.add_job_options({'Variant_calling': {'skip_variants': self._args.skip_variants}})
        if self._args.calling_method:
            self._pipeline.add_job_options({'Variant_calling': {'calling_method': self._args.calling_method}})
        if self._args.output_all_sites:
            self._pipeline.add_job_options({'Variant_calling': {'variants_only': False}})
        if self._args.mutation_rate:
            self._pipeline.add_job_options({'Variant_calling': {'mutation_rate': self._args.mutation_rate}})

    def __save_output_file(self):
        """
        Saves the output file.
        :return: None
        """
        with open(self._pipeline.outputs['VCF'][0].path) as vcf_content:
            with open(self._args.output, 'w') as handle:
                handle.write(vcf_content.read())

if __name__ == '__main__':
    main = SamtoolsVariantCalling()
    main.call_variants()
