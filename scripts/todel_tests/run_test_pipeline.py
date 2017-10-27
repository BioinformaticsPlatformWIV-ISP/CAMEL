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
    DB_LOGGING = False
    SNAKEFILE = os.path.join(os.path.dirname(__file__), 'test.snakefile')

    def __init__(self):
        """
        Initialise.
        :return: None
        """
        self._args = None
        self.config_data = dict()
        self.pipeline = None

        self.runtime_config_name = os.path.join(os.getcwd(),'runtime_config.yaml')


    def __parse_command_line(self):
        """
        Parses the command line arguments.
        :return: Arguments
        """
        args_parser = argparse.ArgumentParser(description='Run the GATK somatic variant caller pipeline.')
        args_parser.add_argument('-W','--wdir', dest='work_dir', metavar='work_dir', help='Working directory', required=True)
        args_parser.add_argument('-R','--fasta_ref', dest='fasta_ref', metavar='fasta_ref', help='Human genome reference fasta file name (as in db_loc).', required=True)
        args_parser.add_argument('-S','--vcf_snps', metavar='vcf_snps', dest='vcf_known_snps', help='Known variant sites (snps) vcf file name (as in db_loc).', required=True)
        args_parser.add_argument('-I','--vcf_indels', metavar='vcf_indels', dest='vcf_known_indels', help='Known variant sites (indels) vcf file name (as in db_loc).', required=True)
        group = args_parser.add_mutually_exclusive_group(required=True)
        group.add_argument('-PE', '--Paired_end', metavar='fq_file', dest='paired_end', help='Paired-end fastq files.', nargs='+')
        group.add_argument('-SE', '--Single_end', metavar='fq_file', dest='single_end', help='Single-end fastq files.', nargs='+')
        args_parser.add_argument('--MarkDuplicates', dest='markduplicates', help='Mark duplicate reads.', action='store_true')

        return args_parser.parse_args()


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
        if self.DB_LOGGING:
            self.pipeline.set_initial_input({'FASTQ_PE': [ToolIOFile(f) for f in self.config_data['fastq_pe']]})

        # Execute the snakemake workflow
        to_execute = 'snakemake --configfile {} --snakefile {}'.format(
            self.runtime_config_name, self.SNAKEFILE
        )
        command = Command(to_execute)
        command.run_command(os.getcwd(), subprocess.STDOUT)
        if command.returncode != 0:
            raise RuntimeError("Error executing Snakemake. Check log for more information")
        # print('Stdout: {}\n'.format(command.stdout))
        # print('Stderr: {}\n'.format(command.stderr))


    def __generate_config_file(self):
        """
        Generates a yaml config file based on CLA.
        :return: None
        """

        print (self._args)

        # Add the job id to the config
        self.config_data['pipeline_job_id'] = self.pipeline.job_id
        self.config_data['pipeline_name'] = self.pipeline.name
        self.config_data['logging'] = self.DB_LOGGING

        self.config_data['working_dir'] = self._args.work_dir
        self.config_data['fasta_ref'] = self._args.fasta_ref
        self.config_data['vcf_known_snps'] = self._args.vcf_known_snps
        self.config_data['vcf_known_indels']=self._args.vcf_known_indels

        if self._args.paired_end:
            self.config_data['fastq_pe'] =self._args.paired_end

        if self._args.single_end:
            self.config_data['fastq_se'] = self._args.single_end

        self.config_data['run_markDuplicates'] = self._args.markduplicates

        # Create config file
        with open(self.runtime_config_name, 'w') as handle:
            yaml.dump(self.config_data, handle)





if __name__ == '__main__':
    main = GATKSomaticMain()
    main.run()