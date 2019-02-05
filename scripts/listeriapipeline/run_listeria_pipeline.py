import argparse
import logging
import os
import sys
import subprocess
import time
import yaml

from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.components.files.fastqutils import FastqUtils
from camel.app.pipeline.pipeline import Pipeline


class ListeriaMain(object):

    """
    Main class to run the Listeria pipeline.
    """

    PIPELINE_VERSION = '1.0'
    SNAKE_FILE = os.path.join(os.path.dirname(__file__), 'pipeline_listeria.snakefile')
    LOGGING_LEVEL = 'pipeline'   # 'step'
    THREADS = 8
    DEBUG = False
    DEBUG_DIR_ROOT = '/scratch/listeria_pipeline/Galaxy_runs/'
    DB_KRAKEN = '/db/kraken/latest/abfhpv'
    # map gene detection db to its db_loc entry (for db path retrieve)
    DB_GENE_DETECTION = {
        'resfinder': 'resfinder',
        'card': 'card',
        'argannot': 'arg_annot',
        'virulencefinder': 'virulencefinder_listeria',
        'plasmidfinder': 'plasmidfinder_grampositive'
    }
    DB_GENE_DETECTION_DB_PARAMS = {
        'resfinder': {
            'min_percent_identity': 90.0,
            'min_coverage': 60.0,
            'version': 'latest',
            'extra_column': ['Antibiotics', 'antibiotics']
        },
        'card': {
            'min_percent_identity': 90.0,
            'min_coverage': 60.0,
            'version': 'latest',
        },
        'argannot': {
            'min_percent_identity': 90.0,
            'min_coverage': 60.0,
            'version': 'latest',
        },
        'virulencefinder': {
            'min_percent_identity': 90.0,
            'min_coverage': 60.0,
            'version': 'latest',
            'extra_column': ['Protein function', 'protein_function']
        },
        'plasmidfinder': {
            'min_percent_identity': 90.0,
            'min_coverage': 60.0,
            'version': 'latest',
            'extra_column': ['Notes', 'notes']
        }
    }
    DB_SEQUENCE_TYPING_ROOT = '/data/sequence_typing/listeria/'
    DB_SEQUENCE_TYPING = {
        'species_confirmation': os.path.join(DB_SEQUENCE_TYPING_ROOT, 'species_confirmation'),
        'MLST-Pasteur': os.path.join(DB_SEQUENCE_TYPING_ROOT, 'mlst'),
        'cgMLST': os.path.join(DB_SEQUENCE_TYPING_ROOT, 'cgmlst'),
        'serogroup': os.path.join(DB_SEQUENCE_TYPING_ROOT, 'serogroup'),
        'virulence': os.path.join(DB_SEQUENCE_TYPING_ROOT, 'virulence'),
        'antibiotic_resistance': os.path.join(DB_SEQUENCE_TYPING_ROOT, 'antibiotic_resistance'),
        'metal_detergent_resistance': os.path.join(DB_SEQUENCE_TYPING_ROOT, 'metal_detergent_resistance')
    }

    def run(self):
        """
        Runs the pipeline.
        :return: None
        """
        self._args = ListeriaMain._parse_arguments()
        self.__preprocess()
        config_file_path = self.__initialize_run_config()
        command = Command('snakemake --cores {} --configfile {} --snakefile {} 2>&1'.format(
            ListeriaMain.THREADS, config_file_path, ListeriaMain.SNAKE_FILE))
        command.run_command(os.getcwd(), subprocess.STDOUT)
        if command.returncode != 0:
            print(command.stdout)
            print(command.stderr)
            raise RuntimeError("Error executing Snakemake. Check log for more information")
        print(command.stdout)
        print(command.stderr)

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
        parser.add_argument('--pubmlst_antibiotic', action='store_true')
        parser.add_argument('--pubmlst_virulence', action='store_true')
        return parser.parse_args()

    def __preprocess(self):
        """
        preprocess based on arguments and prepare to run pipeline
        """
        # sample name
        if self._args.sample_name:
            self._sample_name = self._args.sample_name
        else:
            self._sample_name = FastqUtils.get_sample_name(self._args.fastq_names[0])

        # input data file
        self._fastq_input = self.__symlink_input_files()

        # output directory and files
        if not os.path.isdir(self._args.output_dir):
            os.mkdir(self._args.output_dir)
        if os.path.isfile(self._args.output_html):
            os.remove(self._args.output_html)
        if os.path.isfile(self._args.output_summary):
            os.remove(self._args.output_summary)

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

    def __initialize_run_config(self):
        """
        Initialize run setting and config based on self.DEBUG tag
        :return: config_file_path (string)
        """
        if self.DEBUG:
            logging.info("Run pipeline with debugging settings ...")

            # DEBUG running setting
            self._args.assembler = 'VelvetOptimiser'
            self._args.cgmlst = False
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
            self._args.kraken_db = self.DB_KRAKEN
            self._args.working_dir = os.path.abspath(os.getcwd())
            config_file_path = os.path.join(self._args.working_dir, 'snake_conf.yml')
            self.__create_config_file(config_file_path)

        # Output configuration file
        with open(config_file_path, 'r') as cfg:
            print("-- Running configuration:")
            print(cfg.read())
            print("-- End of configuration")
        sys.stdout.flush()

        return config_file_path

    def __create_config_file(self, path):
        """
        Creates a config file based on the user input.
        :param path: Path
        :return: None
        """
        logging.info("Creating config file.")
        pipeline_name = 'Listeria Pipeline'
        _pipeline = Pipeline(name=pipeline_name, camel=Camel(), logging_level=self.LOGGING_LEVEL)

        with open(path, 'w') as handle:
            yaml.dump({
                'sample_name': self._sample_name,
                'fastq_pe': self._fastq_input,
                'assembler': self._args.assembler,
                'detection_method': self._args.analysis_type,
                'library': self._args.library
            }, handle, default_flow_style=False)
            yaml.dump({
                'pipeline_version': self.PIPELINE_VERSION,
                'pipeline_name': pipeline_name,
                'pipeline_job_id': _pipeline.job_id,
                'logging_level': self.LOGGING_LEVEL,
            }, handle, default_flow_style=False)
            yaml.dump({
                'report_html': self._args.output_html,
                'report_summary': self._args.output_summary,
                'output_dir': self._args.output_dir,
                'working_dir': self._args.working_dir
            }, handle, default_flow_style=False)

            # Gene detection
            yaml.dump({'gene_detection': self.__get_gene_detection_db_config()}, handle, default_flow_style=False)
            # Sequence typing
            yaml.dump({'sequence_typing': self.__get_sequence_typing_db_config()}, handle, default_flow_style=False)

            # Kraken databaxse
            yaml.dump({'db_kraken': self._args.kraken_db}, handle, default_flow_style=False)
            yaml.dump({'kraken_expspecies': 'Listeria monocytogenes'}, handle, default_flow_style=False)

            # SequenceTyping imperfect hits handling: two options,
            # - st_mark_impeefect_hit: imperfect hits is reported with mark '(p)
            # yaml.dump({'st_mark_imperfect_hit': True}, handle, default_flow_style=False)
            # - st_imperfect_as_nohit: imperfect hits will be skip, nohit will be reported in place (pubmlst style)
            yaml.dump({'st_imperfect_as_nohit': True}, handle, default_flow_style=False)

    def __get_gene_detection_db_config(self):
        """
        Returns the gene detection related databases.
        :return: dict of database configs
        """
        dbs = {}
        if self._args.argannot:
            dbs[self.DB_GENE_DETECTION['argannot']] = self.DB_GENE_DETECTION_DB_PARAMS['argannot']
        if self._args.card:
            dbs[self.DB_GENE_DETECTION['card']] = self.DB_GENE_DETECTION_DB_PARAMS['card']
        if self._args.resfinder:
            dbs[self.DB_GENE_DETECTION['resfinder']] = self.DB_GENE_DETECTION_DB_PARAMS['resfinder']
        if self._args.virulencefinder:
            dbs[self.DB_GENE_DETECTION['virulencefinder']] = self.DB_GENE_DETECTION_DB_PARAMS['virulencefinder']
        if self._args.plasmidfinder:
            dbs[self.DB_GENE_DETECTION['plasmidfinder']] = self.DB_GENE_DETECTION_DB_PARAMS['plasmidfinder']
        return dbs

    def __get_sequence_typing_db_config(self):
        """
        Returns the sequence typing related databases
        :return: dictionary of seqtype db and location(path)
        """
        sequence_typing_dbs = {}
        if self._args.species_confirmation:
            sequence_typing_dbs['species_confirmation'] = self.DB_SEQUENCE_TYPING['species_confirmation']
        if self._args.mlst:
            sequence_typing_dbs['MLST-Pasteur'] = self.DB_SEQUENCE_TYPING['MLST-Pasteur']
        if self._args.cgmlst:
            sequence_typing_dbs['cgMLST'] = self.DB_SEQUENCE_TYPING['cgMLST']
        if self._args.serogrouping:
            sequence_typing_dbs['serogroup'] = self.DB_SEQUENCE_TYPING['serogroup']
        if self._args.pubmlst_virulence:
            sequence_typing_dbs['virulence'] = self.DB_SEQUENCE_TYPING['virulence']
        if self._args.pubmlst_antibiotic:
            sequence_typing_dbs['antibiotic_resistance'] = self.DB_SEQUENCE_TYPING['antibiotic_resistance']
        if self._args.pubmlst_metal:
            sequence_typing_dbs['metal_detergent_resistance'] = self.DB_SEQUENCE_TYPING['metal_detergent_resistance']

        return sequence_typing_dbs


if __name__ == '__main__':

    # initialize camel object to initialize logging configuration
    camel = Camel()

    logging.basicConfig(level=logging.DEBUG)
    main = ListeriaMain()
    main.run()
