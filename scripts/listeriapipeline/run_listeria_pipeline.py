import argparse
import logging
import os
import sys
import subprocess
import time
import yaml

from app.command.command import Command
from app.components.files.fastqutils import FastqUtils


class ListeriaMain(object):

    """
    Main class to run the Listeria pipeline.
    """

    PIPELINE_VERSION = '0.1'
    SNAKE_FILE = os.path.join(os.path.dirname(__file__), 'pipeline_listeria.snakefile')
    LOG_IO = False
    THREADS = 8
    DEBUG = False
    DEBUG_DIR_ROOT = '/scratch/qiafu/listeria_pipeline/Galaxy_runs/'

    @staticmethod
    def _parse_arguments():
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        parser.add_argument('--sample-name', help="Name of the sample")
        parser.add_argument('--fastq-pe', nargs=2, help="FASTQ input files")
        parser.add_argument('--fastq-names', nargs=2, help="FASTQ input file names")
        parser.add_argument('--assembler', help="Assembler to use.", choices=['VelvetOptimiser', 'SPAdes'])
        parser.add_argument('--analysis-type', help="Type of analysis (fast / normal)", choices=['fast', 'normal'])
        parser.add_argument('--library', help="Adapter library that was used for the sequencing")
        parser.add_argument('--output-summary', help="Output file for the summary")
        parser.add_argument('--output-html', help="Output file for the HTML report")
        parser.add_argument('--output-dir', help="Output directory for the files in the HTML report")
        # gene detection
        parser.add_argument('--resfinder', action='store_true')
        parser.add_argument('--argannot', action='store_true')
        parser.add_argument('--card', action='store_true')
        parser.add_argument('--virulencefinder', action='store_true')
        parser.add_argument('--plasmidfinder', action='store_true')
        # sequence typing
        parser.add_argument('--species_confirmation', action='store_true')
        parser.add_argument('--mlst', action='store_true')
        parser.add_argument('--cgmlst', action='store_true')
        parser.add_argument('--serogrouping', action='store_true')
        parser.add_argument('--pubmlst_metal', action='store_true')
        parser.add_argument('--pubmlst_resistance', action='store_true')
        parser.add_argument('--pubmlst_virulence', action='store_true')
        return parser.parse_args()

    def __init__(self):
        """
        Initializes the main class.
        """
        self._args = ListeriaMain._parse_arguments()
        if self._args.sample_name:
            self._sample_name = self._args.sample_name
        else:
            self._sample_name = FastqUtils.get_sample_name(self._args.fastq_names[0])
        self._fastq_input = self.__symlink_input_files()

    def __symlink_input_files(self):
        """
        Creates a symbolic link to the input files.
        :return: Complete path to the symlinks
        """
        logging.info("Creating symbolic links for the input files")
        if not os.path.isdir('input_files'):
            os.mkdir('input_files')
        links = []
        for name, file_ in zip(self._args.fastq_names, self._args.fastq_pe):
            path = os.path.abspath(os.path.join('input_files', name))
            os.symlink(file_, path)
            links.append(path)
        return links

    def run(self):
        """
        Runs the pipeline.
        :return: None
        """
        config_file_path = self.__initialize_run_config()

        if not os.path.isdir(self._args.output_dir):
            os.mkdir(self._args.output_dir)
        if os.path.isfile(self._args.output_html):
            os.remove(self._args.output_html)
        command = Command('snakemake --cores {} --configfile {} --snakefile {}'.format(
            ListeriaMain.THREADS, config_file_path, ListeriaMain.SNAKE_FILE))
        command.run_command(os.getcwd(), subprocess.STDOUT)
        if command.returncode != 0:
            print(command.stdout)
            print(command.stderr)
            raise RuntimeError("Error executing Snakemake. Check log for more information")
        print(command.stdout)
        print(command.stderr)

    def __initialize_run_config(self):
        """
        Initialize run setting and config based on self.DEBUG tag
        :return: config_file_path (string)
        """
        if self.DEBUG:
            logging.info("Run pipeline with debugging settings ...")

            # DEBUG running setting
            self._args.assembler = 'VelvetOptimiser'
            self._args.kraken_db = '/data/kraken/latest/abfhpv_lite/'

            # DEBUG report and directory structure
            self._args.working_dir = os.path.join(self.DEBUG_DIR_ROOT, time.strftime("%Y%m%d_%H%M%S"))
            self._args.output_dir = os.path.join(self._args.working_dir, 'final_report')
            os.makedirs(self._args.output_dir)
            # create a softlink to Galaxy report/summary
            sl_report = os.path.join(self._args.output_dir, 'galaxy_run_report.html')
            os.symlink(self._args.output_html, sl_report)
            sl_summary = os.path.join(self._args.output_dir, 'galaxy_summary.tsv')
            os.symlink(self._args.output_summary, sl_summary)
            # Move to working dir, so that everything is created under (ep. for camel.log)
            os.chdir(self._args.working_dir)

            # DEBUG config file
            config_file_path = os.path.join(self._args.working_dir, 'snake_conf.yml')
            self.__create_config_file(config_file_path)

        else:
            self._args.kraken_db = '/data/kraken/latest/abfhpv/'
            self._args.working_dir = os.path.abspath(os.getcwd())
            config_file_path = os.path.join(self._args.working_dir, 'snake_conf.yml')
            self.__create_config_file(config_file_path)

        # Output configuration file
        with open(config_file_path, 'r') as cfg:
            print("-- Running configuration:")
            print(cfg.read())
            print("-- End of  configuration")
        sys.stdout.flush()

        return config_file_path

    def __create_config_file(self, path):
        """
        Creates a config file based on the user input.
        :param path: Path
        :return: None
        """
        logging.info("Creating config file.")
        with open(path, 'w') as handle:
            yaml.dump({'logging': False, 'pipeline_name': 'Listeria Pipeline', 'pipeline_job_id': 3},
                      handle, default_flow_style=False)
            yaml.dump({'report': self._args.output_html}, handle, default_flow_style=False)
            yaml.dump({'output_dir': self._args.output_dir}, handle, default_flow_style=False)
            yaml.dump({'working_dir': self._args.working_dir}, handle, default_flow_style=False)
            yaml.dump({'sample_name': self._sample_name}, handle, default_flow_style=False)
            yaml.dump({'assembler': self._args.assembler}, handle, default_flow_style=False)
            yaml.dump({'detection_method': self._args.analysis_type}, handle, default_flow_style=False)
            yaml.dump({'fastq_pe': self._fastq_input}, handle, default_flow_style=False)

            # Gene detection
            gene_detection_dbs = []
            resistance_dbs = self.__get_resistance_db_config()
            if len(resistance_dbs) > 0:
                gene_detection_dbs = resistance_dbs
            if self._args.virulencefinder:
                gene_detection_dbs.append('VirulenceFinder')
            if self._args.plasmidfinder:
                gene_detection_dbs.append('Gram_positive')
            yaml.dump({'gene_detection': gene_detection_dbs}, handle, default_flow_style=False)

            # Sequence typing
            sequence_typing_dbs = []
            if self._args.species_confirmation:
                sequence_typing_dbs.append('species_confirmation')
            if self._args.mlst:
                sequence_typing_dbs.append('MLST-Pasteur')
            if self._args.cgmlst:
                sequence_typing_dbs.append('cgMLST')
            if self._args.serogrouping:
                sequence_typing_dbs.append('serogroup')
            if self._args.pubmlst_virulence:
                sequence_typing_dbs.append('virulence')
            if self._args.pubmlst_resistance:
                sequence_typing_dbs.append('antibiotic_resistance')
            if self._args.pubmlst_metal:
                sequence_typing_dbs.append('metal_detergent_resistance')
            yaml.dump({'sequence_typing': sequence_typing_dbs}, handle, default_flow_style=False)

            # Kraken databaxse
            yaml.dump({'db_kraken': self._args.kraken_db}, handle, default_flow_style=False)

    def __get_resistance_db_config(self):
        """
        Returns the resistance characterization database config.
        :return: Database config
        """
        dbs = []
        if self._args.argannot:
            dbs.append('ARG-ANNOT')
        if self._args.card:
            dbs.append('CARD')
        if self._args.resfinder:
            dbs.append('ResFinder')
        return dbs


if __name__ == '__main__':

    logging.basicConfig(level=logging.DEBUG)
    main = ListeriaMain()
    main.run()
