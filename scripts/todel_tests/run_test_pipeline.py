import os
import yaml
import argparse

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

    def __init__(self):
        """
        Initialise.
        :return: None
        """
        self.args = self.__parse_command_line()



    def __parse_command_line(self):
        """
        Parses the command line arguments.
        :return: Arguments
        """
        args_parser = argparse.ArgumentParser(description='Run the GATK somatic variant caller pipeline.')
        args_parser.add_argument('-W','--wdir', dest='wdir', metavar='work_dir', help='Working directory', required=True)
        args_parser.add_argument('-R','--fasta_ref', dest='fasta_ref', metavar='fasta_ref', help='Human genome reference fasta file name (as in db_loc).', required=True)
        args_parser.add_argument('-S','--vcf_snps', metavar='vcf_snps', dest='vcf_snps', help='Known variant sites (snps) vcf file name (as in db_loc).', required=True)
        args_parser.add_argument('-I','--vcf_indels', metavar='vcf_indels', dest='vcf_indels', help='Known variant sites (indels) vcf file name (as in db_loc).', required=True)
        group = args_parser.add_mutually_exclusive_group(required=True)
        group.add_argument('-PE', '--Paired_end', metavar='fq_file', dest='paired_end', help='Paired-end fastq files.', nargs='+')
        group.add_argument('-SE', '--Single_end', metavar='fq_file', dest='single_end', help='Single-end fastq files.', nargs='+')
        args_parser.add_argument('--MarkDuplicates', dest='markduplicates', help='Mark duplicate reads.', action='store_true')

        return args_parser.parse_args()


    def run(self):

        # Create a pipeline object
        camel = Camel()
        pipeline = SnakePipeline('GATK somatic calling', camel, self.DB_LOGGING)

        # Load the data of the other config file
        with open(os.path.join(os.path.dirname(__file__), 'config.yaml')) as handle:
            config_data = yaml.load(handle)


        # Add the job id to the config
        config_data['pipeline_job_id'] = pipeline.job_id
        config_data['pipeline_name'] = pipeline.name
        config_data['logging'] = self.DB_LOGGING

        # Setting the initial input makes sure that they are logged
        if self.DB_LOGGING:
            pipeline.set_initial_input({'FASTQ_PE': [ToolIOFile(f) for f in config_data['fastq_pe']]})

        # Create a new config file
        with open('runtime_config.yaml', 'w') as handle:
            yaml.dump(config_data, handle)

        # Execute the snakemake workflow
        command = Command('snakemake --configfile {} --snakefile {}'.format(
            'runtime_config.yaml', os.path.join(os.path.dirname(__file__), 'test.snakefile')
        ))
        command.run_command('.')
        print('Stdout: {}\n'.format(command.stdout))
        print('Stderr: {}\n'.format(command.stderr))



if __name__ == '__main__':
    main = GATKSomaticMain()
    main.run()