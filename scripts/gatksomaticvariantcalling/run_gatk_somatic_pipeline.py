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


class GATKSomaticMain(object):
    """
    Main class to run the GATK somatic variant caller pipeline.
    Generates a config yml file based on CL arguments and runs the pipeline.
    """
    DB_LOGGING = True
    # DEBUG = True
    SNAKEFILE = os.path.join(os.path.dirname(__file__), 'gatk_somatic_steps.snakefile')
    # FROM_GALAXY = False
    CORES = 5

    def __init__(self):
        """
        Initialise.
        :return: None
        """
        self._args = None
        self._config_data = dict()
        self._pipeline = None

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

        # Variant caller to use
        ap.add_argument('-V', '--variant_caller', metavar='variant_caller', dest='variant_caller', help='Variant caller to use.', choices=["mutect1", "mutect2"])

        # output
        ap.add_argument('--mutect1_vcf_output', dest='mutect1_vcf_output', metavar='mutect1_vcf_output', help='Output vcf file from MuTect1.')
        ap.add_argument('--mutect1_tab_output', dest='mutect1_tab_output', metavar='mutect1_tab_output', help='Output variant-call table file for MuTect1.')
        ap.add_argument('--mutect2_vcf_output', dest='mutect2_vcf_output', metavar='mutect2_vcf_output', help='Output vcf file from MuTect2.')

        ap.add_argument('--covar_output', dest='covar_output', metavar='covar_output', help='Output covariates analysis pdf')
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
        ap.add_argument('--downsampling_type', dest='downsampling_type',
                        help='Type of downsampling to performe on reads (by MuTect/MuTect2). NONE,ALL_READS,BY_SAMPLE. Default: BY_SAMPLE. Perform or not downsampling on reads. '
                             'By default, MuTect/MuTect2 downsamples to 1000 reads. Usage example: --downsample None (disables downsampling).')
        ap.add_argument('--downsampling_target', dest='downsampling_target', help='Target value for downsampling to perform on reads (by MuTect/MuTect2). Default: 1000. Usage example: --downsample 10000 (sets target value to 10000 reads).')
        # MuTect1:
        # gap_events_threshold
        ap.add_argument('--gap_events_threshold', dest='gap_events_threshold', help='For MuTect1; number of reads allowed to contain insdels around a fixed window (MuTect1 default 11 bp) before being marked as gap_event and filtered-out.')

        # strand_artifact_lod
        ap.add_argument('--strand_artifact_lod', dest='strand_artifact_lod', help='For MuTect1; log-odds ratio for strand bias. Default MuTect: 2.0; disable: -99999')

        # MuTect2:
        # run from galaxy flag
        ap.add_argument('--from_galaxy', dest='from_galaxy', help='Indicates that the command is run from galaxy. Useful for logging stderr.', action='store_true')

        # number of threads
        ap.add_argument('--threads', dest='threads', help='Maximum number of threads for snakemake to allow.', default=self.CORES)

        # snakemake arguments
        # snakemake unlock
        ap.add_argument('--unlock', dest='unlock', action="store_const", const="--unlock", help='Unlock snakemake working directory.', default="")
        # snakemake dag
        ap.add_argument('--dag', dest='dag', action="store_const", const="--dag", help='Generate snakemake DAG.', default="")
        # dryrun
        # snakemake dag
        ap.add_argument('--dryrun', dest='dryrun', action="store_const", const="--dryrun", help='Snakemake dryrun; generates the list of jobs.', default="")

        # job id
        ap.add_argument('--job_id', dest='job_id', metavar='job_id', help='Job ID for debugging and logging.',
                        default=datetime.datetime.now().strftime("%Y%m%d_%H%M%S-%f"))

        return ap.parse_args()

    def __generate_config_file(self):
        """
        Generates a yaml config file based on CLA.
        :return: None
        """

        # Name of config file generated at runtime for snakemake pipeline
        self.runtime_config_name = os.path.join(os.getcwd(), 'runtime_config_{}.yaml'.format(self._args.job_id))

        # Add the job id to the config
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

        # variant caller choice
        self._config_data['variant_caller'] = self._args.variant_caller

        # Output filenames
        if self._args.mutect1_vcf_output:
            self._config_data['mutect1_vcf_output'] = self._args.mutect1_vcf_output
        if self._args.mutect1_tab_output:
            self._config_data['mutect1_tab_output'] = self._args.mutect1_tab_output
        if self._args.mutect2_vcf_output:
            self._config_data['mutect2_vcf_output'] = self._args.mutect2_vcf_output
        if self._args.covar_output:
            self._config_data['covar_output'] = self._args.covar_output
        if self._args.bam_output:
            self._config_data['bam_output'] = self._args.bam_output

        # Flag for MarkDuplicates
        self._config_data['run_markDuplicates'] = self._args.markduplicates

        # MuTect1 parameters
        if self._config_data['variant_caller'] == "mutect1":
            # Downsampling
            if self._args.downsampling_type:
                self._config_data['MuTect1_downsampling_type'] = self._args.downsampling_type
            if self._args.downsampling_target:
                self._config_data['MuTect1_downsampling_target'] = self._args.downsampling_target
            # gap_event_threshold
            if self._args.gap_events_threshold:
                self._config_data['gap_events_threshold'] = self._args.gap_events_threshold
            # strand_artifact_lod
            if self._args.strand_artifact_lod:
                self._config_data['strand_artifact_lod'] = self._args.strand_artifact_lod

        # MuTect2 parameters
        if self._config_data['variant_caller'] == "mutect2":
            # Downsampling
            if self._args.downsampling_target:
                self._config_data['MuTect2_downsampling_target'] = self._args.downsampling_target

        # Indel realigner flag (indel realignment only required for MuTect1)
        if self._args.variant_caller == "mutect1":
            self._config_data['run_indel_realignment'] = True
        elif self._args.variant_caller == "mutect2":
            self._config_data['run_indel_realignment'] = False

        # threads
        if self._args.threads:
            self._config_data['threads'] = self._args.threads
        else:
            self._config_data['threads'] = self._args.threads
        # Create and write to config file
        with open(self.runtime_config_name, 'w') as handle:
            yaml.dump(self._config_data, handle)

    def run(self):
        """
        Sets-up and runs the pipeline and logs stderr and stdout if running the command fails.
        :return: None 
        """

        # Create a pipeline object
        camel = Camel()
        self._pipeline = SnakePipeline('GATK somatic calling', camel, self.DB_LOGGING)

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
        snakemake_params = " ".join([self._args.unlock, self._args.dag, self._args.dryrun])
        to_execute = 'snakemake --configfile {config} --snakefile {snakefile} --cores {cores} {snakemake_params}'.format(config=self.runtime_config_name, snakefile=self.SNAKEFILE, cores=self._args.threads, snakemake_params=snakemake_params)
        command = Command(to_execute)
        command.run_command(self._args.work_dir, subprocess.STDOUT)
        if command.returncode != 0 and self._args.from_galaxy:
            with open("/scratch/temp/galaxy_logs/{}_Stdout".format(self._args.job_id), "w") as file_out:
                file_out.write(command.stdout)
            with open("/scratch/temp/galaxy_logs/{}_Stderr".format(self._args.job_id), "w") as file_out:
                file_out.write(command.stderr)
            raise RuntimeError(
                "Error executing Snakemake. Check log ('/scratch/temp/galaxy_logs/{}_Stderr') for more information.".format(
                    self._args.job_id))


if __name__ == '__main__':
    main = GATKSomaticMain()
    main.run()
