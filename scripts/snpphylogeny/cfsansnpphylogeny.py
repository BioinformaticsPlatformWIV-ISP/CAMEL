import os

from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue
from app.tools.cfsan.cfsansnppipeline import CfsanSnpPipeline
from scripts.snpphylogeny.cfsan import CFSAN_TEMPLATE
from scripts.snpphylogeny.snpphylogeny import SnpPhylogeny


class CFSANSnpPhylogeny(SnpPhylogeny):
    """
    SNP phylogeny pipeline based on the CFSAN SNP pipeline.
    """

    BOWTIE2_THREADS = 8
    CONCURRENT_SAMPLES = 4

    def __init__(self):
        """
        Initializes the pipeline script.
        """
        super(CFSANSnpPhylogeny, self).__init__('CFSAN')

    @staticmethod
    def _add_specific_arguments(argument_parser):
        """
        Adds the pipeline specific arguments.
        :param argument_parser: Argument parser
        :return: None
        """
        argument_parser.add_argument('--selected-matrix', choices=['preserved', 'regular'], required=True)

    def _run_snp_calling(self, reads):
        """
        Runs the SNP calling on the input samples.
        :param reads: List of PE reads for all samples
        :return: SNP matrix, output files
        """
        cfsan = CfsanSnpPipeline(self._camel)

        cfsan.add_input_files({'FASTA': [ToolIOFile(self._args.reference)],
                               'TXT': [ToolIOFile(self.__create_config_file())]})
        for input_set, sample_name in zip(reads, self._sample_names):
            cfsan.add_input_files({'FASTQ': input_set['FASTQ_PE']})
            cfsan.add_input_files({'VAL_Name': [ToolIOValue(sample_name)]})
        cfsan.run(os.path.abspath('.'))
        self._cfsan_informs = cfsan.informs

        if self._args.selected_matrix == 'preserved':
            snp_matrix = cfsan.tool_outputs['FASTA_Preserved'][0]
        else:
            snp_matrix = cfsan.tool_outputs['FASTA'][0]
        return snp_matrix, CFSANSnpPhylogeny.__get_output_files(cfsan.tool_outputs)

    @staticmethod
    def __create_config_file():
        """
        Creates the config file based on the template.
        :return: Path to config file
        """
        with open(CFSAN_TEMPLATE) as input_handle:
            template_data = input_handle.read()
        config_path = os.path.join(os.path.abspath('.'), 'cfsan.conf')
        with open(config_path, 'w') as output_handle:
            output_handle.write(template_data.format(bowtie2_threads=CFSANSnpPhylogeny.BOWTIE2_THREADS,
                                                     concurrent_samples=CFSANSnpPhylogeny.CONCURRENT_SAMPLES))
        return config_path

    @staticmethod
    def __get_output_files(cfsan_outputs):
        """
        Returns the output files
        :param cfsan_outputs: CFSAN outputs
        :return: List with: (Name, Filename template, Output files)
        """
        return [
            ('Alignment (BAM)', 'alignment_{}.bam', cfsan_outputs['BAM']),
            ('SNPs filtered (VCF)', 'snps_{}.vcf', cfsan_outputs['VCF']),
            ('SNPs filtered preserved (VCF)', 'snps_preserved_{}.vcf', cfsan_outputs['VCF_preserved'])
        ]

    COLUMN_MAPPING = [
        ('Average_Insert_Size', 'Insert Size (avg)'),
        ('Average_Pileup_Depth', 'Pileup Depth (avg)'),
        ('Percent_of_Reads_Mapped', '% Mapped'),
        ('Phase1_Preserved_SNPs', 'SNPs pres. (P1)'),
        ('Phase1_SNPs', 'SNPs (P1)'),
        ('Phase2_Preserved_SNPs', 'SNPs pres. (P2)'),
        ('Phase2_SNPs', 'SNPs (P2)')
    ]

    def _collect_metrics(self):
        """
        Collects the per sample analysis metrics.
        :return: List of header + values for each sample
        """
        header = [mapping[1] for mapping in CFSANSnpPhylogeny.COLUMN_MAPPING]
        data = []
        for sample_name in sorted(self._cfsan_informs.keys()):
            data.append([self._cfsan_informs[sample_name][key] for key in [mapping[0] for mapping in
                                                                           CFSANSnpPhylogeny.COLUMN_MAPPING]])
        return [header] + data


if __name__ == '__main__':
    cfsan_snp = CFSANSnpPhylogeny()
    cfsan_snp.run()
