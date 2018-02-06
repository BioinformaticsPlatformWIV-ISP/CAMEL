import os
import yaml
import argparse
import subprocess

from app.camel import Camel
from app.command.command import Command
from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue
from app.pipeline.snakepipeline import SnakePipeline
import datetime
from config import MAIN_CONFIG


class GATKSomaticMain(object):
    """
    Main class to run the GATK somatic variant caller pipeline.
    Generates a config yml file based on CL arguments and runs the pipeline.
    """
    DB_LOGGING = True
    SNAKEFILE = os.path.join(os.path.dirname(__file__), 'gatk_somatic_steps.snakefile')
    CORES = 5

    def __init__(self):
        """
        Initialise.
        :return: None
        """
        self._args = None
        self._config_data = dict()
        self._pipeline = None

        # Name of config file generated at runtime for snakemake pipeline
        self.runtime_config_name = os.path.join(os.getcwd(), 'runtime_config.yaml')

        self.camel = Camel()

        # set galaxy dump directory in case of failure
        self._galaxy_dump_dir = os.path.join(self.camel.config["galaxy"]["dump_dir"], "GATK_somatic_calling")

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
        gp.add_argument('-PE', '--paired_end', metavar='fq_file', dest='paired_end', help='Paired-end fastq files.',
                        nargs='+')
        gp.add_argument('-SE', '--single_end', metavar='fq_file', dest='single_end', help='Single-end fastq files.',
                        nargs='+')

        # output
        ap.add_argument('--vcf_output', dest='vcf_output', metavar='vcf_output', help='Output vcf file')
        ap.add_argument('--tab_output', dest='tab_output', metavar='tab_output', help='Output variant-call table file')
        ap.add_argument('--covar_output', dest='covar_output', metavar='covar_output',
                        help='Output covariates analysis pdf')
        ap.add_argument('--bam_output', dest='bam_output', metavar='bam_output', help='Aligned BAM file')

        # references
        ap.add_argument('-R', '--fasta_ref', dest='fasta_ref', metavar='fasta_ref',
                        help='Human genome reference fasta file name (as in db_loc).', required=True)
        ap.add_argument('-S', '--vcf_snps', metavar='vcf_snps', dest='vcf_known_snps',
                        help='Known variant sites (snps) vcf file name (as in db_loc).', required=True)
        ap.add_argument('-I', '--vcf_indels', metavar='vcf_indels', dest='vcf_known_indels',
                        help='Known variant sites (indels) vcf file name (as in db_loc).', required=True)

        # MarkDuplicates flag
        ap.add_argument('--mark_duplicates', dest='markduplicates', help='Mark duplicate reads.', action='store_true')

        # Downsampling
        ap.add_argument('--downsampling_type', dest='downsampling_type', help='Type of downsampling to performe on reads (by MuTect). NONE,ALL_READS,BY_SAMPLE. Default: BY_SAMPLE. Perform or not downsampling on reads. By default, MuTect downsamples to 1000 reads. Usage example: --downsample None (disables downsampling).')
        ap.add_argument('--downsampling_target', dest='downsampling_target', help='Target value for downsampling to performe on reads (by MuTect). Default: 1000. Usage example: --downsample 10000 (sets target value to 10000 reads).')

        # gap_events_threshold (MuTect)
        ap.add_argument('--gap_events_threshold', dest='gap_events_threshold', help='Number of reads allowed to contain insdels around a fixed window (MuTect default 11 bp) before being marked as gap_event and filtered-out.')

        # gap_events_threshold (MuTect)
        ap.add_argument('--strand_artifact_lod', dest='strand_artifact_lod', help='Log-odds ratio for strand bias. Default MuTect: 2.0; disable: -99999')

        # run from galaxy flag
        ap.add_argument('--from_galaxy', dest='from_galaxy', help='Indicates that the command is run from galaxy. Useful for logging stderr.', action='store_true')

        # job id
        ap.add_argument('--job_id', dest='job_id', metavar='job_id', help='Job ID for debugging and logging. (Not the same as internal camel pipeline job id!)',
                        default=datetime.datetime.now().strftime("%Y%m%d_%H%M%S-%f"))

        return ap.parse_args()

    def __generate_config_file(self):
        """
        Generates a yaml config file based on CLA.
        :return: None
        """

        # Add the camel pipeline info to the config
        self._config_data['pipeline_job_id'] = self._pipeline.job_id
        self._config_data['pipeline_name'] = self._pipeline.name
        self._config_data['logging'] = self.DB_LOGGING

        # Set working directory: if not given as CLA, set to current working directory (Galaxy)
        if self._args.work_dir:
            self._config_data['working_dir'] = self._args.work_dir
        else:
            self._config_data['working_dir'] = os.getcwd()

        # References
        self._config_data['fasta_ref'] = self._args.fasta_ref
        self._config_data['vcf_known_snps'] = self._args.vcf_known_snps
        self._config_data['vcf_known_indels'] = self._args.vcf_known_indels

        # Input fastq files
        if self._args.paired_end:
            self._config_data['fastq'] = self._args.paired_end
            self._config_data['PE'] = True
            self._config_data['SE'] = False
        if self._args.single_end:
            self._config_data['fastq'] = self._args.single_end
            self._config_data['SE'] = True
            self._config_data['PE'] = False

        # Output filenames
        if self._args.vcf_output:
            self._config_data['vcf_output'] = self._args.vcf_output
        if self._args.tab_output:
            self._config_data['txt_output'] = self._args.tab_output
        if self._args.covar_output:
            self._config_data['covar_output'] = self._args.covar_output
        if self._args.bam_output:
            self._config_data['bam_output'] = self._args.bam_output

        # Flag for MarkDuplicates
        self._config_data['run_markDuplicates'] = self._args.markduplicates

        # MuTect parameters
        # Downsampling
        if self._args.downsampling_type:
            self._config_data['downsampling_type'] = self._args.downsampling_type
        if self._args.downsampling_target:
            self._config_data['downsampling_target'] = self._args.downsampling_target
        # gap_event_threshold
        if self._args.gap_events_threshold:
            self._config_data['gap_events_threshold'] = self._args.gap_events_threshold
        # strand_artifact_lod
        if self._args.strand_artifact_lod:
            self._config_data['strand_artifact_lod'] = self._args.strand_artifact_lod

        # Create and write to config file
        with open(self.runtime_config_name, 'w') as handle:
            yaml.dump(self._config_data, handle)

    def run(self):
        """
        Sets-up and runs the pipeline and logs stderr and stdout if running the command fails.
        :return: None 
        """

        # Create a pipeline object

        self._pipeline = SnakePipeline('GATK somatic calling', self.camel, self.DB_LOGGING)

        self._args = self.__parse_command_line()

        self.__generate_config_file()

        # Set the initial input
        input_dict = dict()
        if self.DB_LOGGING:
            if self._args.paired_end:
                input_dict['FASTQ_PE'] = [ToolIOFile(f) for f in self._config_data['fastq']]
            elif self._args.single_end:
                input_dict['FASTQ_SE'] = [ToolIOFile(f) for f in self._config_data['fastq']]
            input_dict['MarkDuplicates'] = [ToolIOValue(self._args.markduplicates)]

            self._pipeline.set_initial_input(input_dict)

        # Execute the snakemake workflow and log stdout and stderr if command fails (if pipeline is run from galaxy).
        to_execute = 'snakemake --configfile {} --snakefile {} --cores {}'.format(self.runtime_config_name, self.SNAKEFILE, self.CORES)
        command = Command(to_execute)
        command.run_command(self._args.work_dir, subprocess.STDOUT)
        if command.returncode != 0 and self._args.from_galaxy:
            with open(os.path.join(self._galaxy_dump_dir, "{}_Stdout".format(self._args.job_id)), "w") as file_out:
                file_out.write(command.stdout)
            with open(os.path.join(self._galaxy_dump_dir, "{}_Stderr".format(self._args.job_id)), "w") as file_out:
                file_out.write(command.stderr)
            raise RuntimeError(
                "Error executing Snakemake. Check log ('{}') for more information.".format(os.path.join(self._galaxy_dump_dir, "{}_Stderr".format(
                    self._args.job_id))))


if __name__ == '__main__':
    main = GATKSomaticMain()
    main.run()
