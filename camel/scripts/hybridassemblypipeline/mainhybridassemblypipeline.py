#!/usr/bin/env python
import argparse
import logging
from pathlib import Path
from typing import Optional, Sequence, List, Dict, Any, Tuple

import pandas as pd
import pkg_resources
import yaml

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.galaxy.galaxyutils import GalaxyUtils
from camel.app.components.pipelines.basepipeline import BasePipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.hybridassemblypipeline import CONFIG_DATA


class MainHybridAssemblyPipeline(BasePipeline):
    """
    Main class to run the hybrid assembly pipeline.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the hybrid assembly pipeline class.
        :param args: arguments to be parsed
        :return: None
        """
        self._data_models = pd.read_table(pkg_resources.resource_filename(
            'camel', 'scripts/hybridassemblypipeline/basecalling_models.tsv'))
        self._args = MainHybridAssemblyPipeline._parse_arguments(args)
        _path_snakefile = pkg_resources.resource_filename('camel', 'scripts/hybridassemblypipeline/snakefile/main.smk')
        super().__init__(
            'Hybrid assembly pipeline', '0.1', _path_snakefile, args)
        if self._args.output_dir is None:
            self._args.output_dir = self._args.output_html.parent
        self._sample_name = MainHybridAssemblyPipeline._determine_sample_name(self._args)

    @staticmethod
    def _determine_sample_name(args: argparse.Namespace) -> str:
        """
        Determines the sample name from the provided arguments.
        :param args: Command line arguments
        :return: Parse sample name
        """
        if args.sample_name is not None:
            return args.sample_name
        logging.debug('Sample name not provided, determining name from short-reads')
        if args.fastq_pe_names is not None:
            return GalaxyUtils.determine_sample_name_from_fq([Path(fq) for fq in args.fastq_pe_names])
        return GalaxyUtils.determine_sample_name_from_fq([Path(fq) for fq in args.fastq_pe])

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()

        # Input
        argument_parser.add_argument('--sample-name', type=str)
        argument_parser.add_argument('--fastq-pe', type=Path, help='Input Fastq PE file', nargs=2, required=True)
        argument_parser.add_argument('--fastq-pe-names', help='Input Fastq PE file', nargs=2)
        argument_parser.add_argument('--fastq-se', type=Path, help='Input Fastq SE files')
        argument_parser.add_argument('--fastq-se-name', help='Input Fastq SE file names')

        # Output
        argument_parser.add_argument('--working-dir', type=Path, default=Path.cwd())
        argument_parser.add_argument('--output-html', type=Path, required=True)
        argument_parser.add_argument('--output-dir', type=Path)

        # Parameters
        argument_parser.add_argument('--ploidy', type=int, choices=[1, 2], default=1)
        argument_parser.add_argument('--ont-qual', type=str, required=True,
                                     choices=['nano-corr', 'nano-hq', 'nano-raw'], default='nano-corr')
        argument_parser.add_argument('--meta', action='store_true', default=None)
        argument_parser.add_argument('--expected-genome-size', type=str, required=True)
        argument_parser.add_argument('--ont-basecalling-model', type=str, choices=[
            'r1041_e82_260bps_fast_g632', 'r1041_e82_260bps_hac_g632', 'r1041_e82_260bps_sup_g632',
            'r1041_e82_400bps_fast_g615', 'r1041_e82_400bps_fast_g632', 'r1041_e82_400bps_hac_g615',
            'r1041_e82_400bps_hac_g632',  'r1041_e82_400bps_sup_g615', 'r104_e81_hac_g5015', 'r104_e81_sup_g5015',
            'r104_e81_sup_g610', 'r941_min_hac_g507', 'r941_min_high_g360', 'r941_min_sup_g507', 'r941_prom_hac_g507',
            'r941_prom_high_g360', 'r941_prom_sup_g507'], default='r941_prom_sup_g507')
        argument_parser.add_argument('--filtlong-keep-percent', type=int)

        # Variant calling
        argument_parser.add_argument('--freebayes-min-alternate-fraction', type=float, default=0.5)
        argument_parser.add_argument('--freebayes-min-alternate-count', type=int, default=10)
        argument_parser.add_argument('--clair3-long-indel', action='store_true', default=None)
        argument_parser.add_argument('--sniffles-min-svlen', type=int)
        argument_parser.add_argument('--sniffles-mapq', type=int)
        argument_parser.add_argument('--sniffles-min-support', type=int)

        # Logging & resources
        argument_parser.add_argument(
            '--galaxy-job-id', type=str, help='Job id of the run in galaxy (used for logging')
        argument_parser.add_argument(
            '--log', action='store_true', help="If this flag is set, config file and error logs are kept")
        argument_parser.add_argument('--threads', type=int, default=8)
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the hybrid assembly pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        path_config = self.__create_snakemake_config_data(input_files)
        self._run_snakemake_main(path_config)

    def _get_fastq_input_links(self) -> List[List[Tuple[Path, str]]]:
        """
        Returns the links to the input FASTQ files.
        :return: Links
        """
        links = []
        input_fastq = self._args.fastq_pe
        input_fastq.extend([self._args.fastq_se])
        for read_nb, path in enumerate(input_fastq, start=1):
            gzipped = FileSystemHelper.is_gzipped(path)
            if read_nb > 2:
                read_nb = 'ont'
            links.append([path, f"{self._sample_name}_{read_nb}.fastq{'.gz' if gzipped else ''}"])
        return links

    def _symlink_input(self) -> List[Dict[str, Any]]:
        """
        Symlinks the input files.
        :return: List of FASTQ input dictionaries
        """
        # Determine link names
        links = self._get_fastq_input_links()

        # Create directory
        dir_links = self._args.working_dir / 'input'
        if not dir_links.exists():
            dir_links.mkdir(parents=True)

        # Link files
        paths_new = []
        for path_orig, link_name in links:
            path_new = dir_links / link_name
            logging.debug(f"Symlinking input file: {path_orig} -> {link_name}")
            if path_new.is_symlink():
                path_new.unlink()
            path_new.symlink_to(path_orig)
            paths_new.append(path_new)

        # Return output dictionary
        return [{'name': p.name, 'path': p} for p in paths_new]

    def __create_snakemake_config_data(self, input_fastq) -> str:
        """
        Creates a Snakemake configuration file.
        :return: Config file data
        """
        config_data = {
            # Pipeline
            'pipeline': {
                'name': self.name,
                'version': self.version,
                'title': 'Hybrid assembly pipeline'
            },
            # Input
            'sample_name': self._sample_name, 'working_dir': self._args.working_dir,
            'input': {
                'illumina': [str(input_fastq[0]['path'].absolute()), str(input_fastq[1]['path'].absolute())],
                'ont': str(input_fastq[2]['path'].absolute())
            },
            # Output
            'output_html': str(self._args.output_html.absolute()),
            'output_dir': str(self._args.output_dir.absolute()),
            # Parameters
            'filtlong': {
                'keep_percent': self._args.filtlong_keep_percent if self._args.filtlong_keep_percent is not None else 95
            },
            'assembly': {
                'flye': {
                    'genome_size': self._args.expected_genome_size,
                    'meta': True if self._args.meta else False,
                    'nano_corr': True if self._args.ont_qual == 'nano-corr' else False,
                    'nano_hq': True if self._args.ont_qual == 'nano-hq' else False,
                    'nano_raw': True if self._args.ont_qual == 'nano-raw' else False},
                'min_contig_length': 1000
            },
            'polishing': {
                'medaka': {
                    'consensus': {'model': self._args.ont_basecalling_model},
                    'stitch': {}
                },
                'polca': {},
                'polypolish': {}
            },
            'read_mapping': {
                'bwa': {},
                'minimap2': {}
            },
            'freebayes': {
                'ploidy': self._args.ploidy,
                'min_alternate_fraction': self._args.freebayes_min_alternate_fraction,
                'min_alternate_count': self._args.freebayes_min_alternate_count
            },
            'clair3': {
                'haploid_precise': True if self._args.ploidy == 1 else False,
                'long_indel': True if self._args.clair3_long_indel is not None else False,
                'model_path': str(self._get_clair3_model(self._args.ont_basecalling_model))
            },
            'sniffles': {
                'mapq': self._args.sniffles_mapq if self._args.sniffles_mapq is not None else 25,
                'min_support':
                    self._args.sniffles_min_support if self._args.sniffles_min_support is not None else 'auto',
                'min_svlen': self._args.sniffles_min_svlen if self._args.sniffles_min_svlen is not None else 35
            }
        }
        with open(CONFIG_DATA) as handle_in:
            mainscriptutils.dict_merge(
                config_data, yaml.safe_load(handle_in.read()))
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    def _get_clair3_model(self, medaka_model: str) -> Path:
        """
        Gets the most suitable clair3 model for the given medaka model.
        From documentation: the medaka model with the highest version equal to or less than the guppy version should
        be selected.
        :param medaka_model: Medaka basecalling model
        :return: Path to clair3 model
        """
        clair3_model = next(
            r['clair3_model'] for r in self._data_models.to_dict('records') if r['medaka_model'] == medaka_model)
        logging.info(f"Best matching Clair3 model for medaka model '{medaka_model}': {clair3_model}")
        path_model = Path(Camel.get_instance().config['db_root'], 'clair3', 'models', clair3_model)
        if not path_model.exists():
            raise ValueError(f'Clair3 model not found: {path_model}')
        return path_model


if __name__ == '__main__':
    Camel.get_instance()
    main = MainHybridAssemblyPipeline()
    main.run()
