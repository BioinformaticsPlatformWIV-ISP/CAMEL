import os
import yaml
import argparse
import subprocess

from app.camel import Camel
from app.command.command import Command
from app.io.tooliofile import ToolIOFile
from app.pipeline.snakepipeline import SnakePipeline



class GATKSomaticMain(object):
    """
    Main class to run the GATK somatic variant caller pipeline.
    Generates a config yml file based on CL arguments and runs the pipeline.
    """
    DB_LOGGING = True
    SNAKEFILE = os.path.join(os.path.dirname(__file__), 'gatk_somatic_steps.snakefile')


    def __init__(self):
        """
        Initialise.
        :return: None
        """
        self._args = None
        self.config_data = dict()
        self.pipeline = None

        # Name of config file generated at runtime for snakemake pipeline
        self.runtime_config_name = os.path.join(os.getcwd(),'runtime_config.yaml')


    def __parse_command_line(self):
        """
        Parses the command line arguments.
        :return: Arguments
        """
        ap = argparse.ArgumentParser(description='Run the GATK somatic variant caller pipeline.')

        # various
        ap.add_argument('-W', '--wdir', dest='work_dir', metavar='work_dir', help='Working directory')

        # input
        gp = ap.add_mutually_exclusive_group(required=True)
        gp.add_argument('-PE', '--Paired_end', metavar='fq_file', dest='paired_end', help='Paired-end fastq files.', nargs='+')
        gp.add_argument('-SE', '--Single_end', metavar='fq_file', dest='single_end', help='Single-end fastq files.', nargs='+')

        # output
        ap.add_argument('--vcf_output', dest='vcf_output', metavar='vcf_output', help='Output vcf file')
        ap.add_argument('--tab_output', dest='tab_output', metavar='tab_output', help='Output variant-call table file')
        ap.add_argument('--covar_output', dest='covar_output', metavar='covar_output', help='Output covariates analysis pdf')
        ap.add_argument('--bam_output', dest='bam_output', metavar='bam_output', help='Aligned BAM file')

        # references
        ap.add_argument('-R','--fasta_ref', dest='fasta_ref', metavar='fasta_ref', help='Human genome reference fasta file name (as in db_loc).', required=True)
        ap.add_argument('-S','--vcf_snps', metavar='vcf_snps', dest='vcf_known_snps', help='Known variant sites (snps) vcf file name (as in db_loc).', required=True)
        ap.add_argument('-I','--vcf_indels', metavar='vcf_indels', dest='vcf_known_indels', help='Known variant sites (indels) vcf file name (as in db_loc).', required=True)

        # MarkDuplicates flag
        ap.add_argument('--mark_duplicates', dest='markduplicates', help='Mark duplicate reads.', action='store_true')

        return ap.parse_args()


    def __generate_config_file(self):
        """
        Generates a yaml config file based on CLA.
        :return: None
        """

        # Add the job id to the config
        self.config_data['pipeline_job_id'] = self.pipeline.job_id
        self.config_data['pipeline_name'] = self.pipeline.name
        self.config_data['logging'] = self.DB_LOGGING

        # Set working directory: if not given as CLA, set to current working directory (Galaxy)
        if not self._args.work_dir:
            self.config_data['working_dir'] = os.getcwd()
        else:
            self.config_data['working_dir'] = self._args.work_dir

        # References
        self.config_data['fasta_ref'] = self._args.fasta_ref
        self.config_data['vcf_known_snps'] = self._args.vcf_known_snps
        self.config_data['vcf_known_indels']=self._args.vcf_known_indels

        # Input fastq files
        if self._args.paired_end:
            self.config_data['fastq'] =self._args.paired_end
            self.config_data['PE'] = True
            self.config_data['SE'] = False
        if self._args.single_end:
            self.config_data['fastq'] = self._args.single_end
            self.config_data['SE'] = True
            self.config_data['PE'] = False

        # Output filenames
        if self._args.vcf_output:
            self.config_data['vcf_output'] = self._args.vcf_output
        if self._args.tab_output:
            self.config_data['txt_output'] = self._args.tab_output
        if self._args.covar_output:
            self.config_data['covar_output'] = self._args.covar_output
        if self._args.bam_output:
            self.config_data['bam_output'] = self._args.bam_output

        # Flag for MarkDuplicates
        self.config_data['run_markDuplicates'] = self._args.markduplicates

        # Create and write to config file
        with open(self.runtime_config_name, 'w') as handle:
            yaml.dump(self.config_data, handle)


    def run(self):
        """
        Sets-up and runs the pipeline.  
        :return: None 
        """

        # Create a pipeline object
        camel = Camel()
        self.pipeline = SnakePipeline('GATK somatic calling', camel, self.DB_LOGGING)

        self._args = self.__parse_command_line()

        self.__generate_config_file()

        # Setting the initial input makes sure that they are logged
        # if self.DB_LOGGING:
        #     self.pipeline.set_initial_input({'FASTQ_PE': [ToolIOFile(f) for f in self.config_data['fastq_pe']]})

        # Execute the snakemake workflow
        to_execute = 'snakemake --configfile {} --snakefile {}'.format(self.runtime_config_name, self.SNAKEFILE)
        command = Command(to_execute)
        command.run_command(os.getcwd(), subprocess.STDOUT)
        if command.returncode != 0:
            print('Stdout: {}\n'.format(command.stdout))
            print('Stderr: {}\n'.format(command.stderr))
            raise RuntimeError("Error executing Snakemake. Check log for more information")


if __name__ == '__main__':
    main = GATKSomaticMain()
    main.run()