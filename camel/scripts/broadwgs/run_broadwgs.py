#!/usr/bin/env python
import os
import argparse
import yaml
import logging
import glob

from pathlib import Path
from typing import Any, Optional, List, Dict, Sequence

from camel.app.camel import Camel
from camel.app.pipeline.pipeline import Pipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.broadwgs import SNAKEFILE_MAIN, CONFIG_DATA


class MainBroadWGSPipeline(object):
    """
    Main class to run the GATK best practices workflow
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        self._name = 'Broad institute WGS GATK Best Practices'
        self._version = '0.5'
        self._snakefile = SNAKEFILE_MAIN
        self._args = self._parse_arguments(args)
        self._working_dir = Path(self._args.working_dir)
        self._pipeline = Pipeline(self._name, Camel.get_instance(), self._args.log, self._args.log)

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)

    def _symlink_input(self) -> List[Dict[str, Any]]:
        """
        Symlinks the input files.
        :return: List of FASTQ input dictionaries
        """
        # Determine link names
        links = glob.glob(str(self._args.ubam_prefix) + "*.unmapped.bam")

        # Create directory
        dir_links = self._working_dir / 'input' / self._args.sample
        if not dir_links.exists():
            dir_links.mkdir(parents=True)

        # Link files
        paths_new = []
        for path_orig in links:
            link_name = os.path.basename(path_orig)
            path_new = os.path.join(dir_links, link_name)
            logging.debug(f"Symlinking input file: {link_name}")
            if os.path.islink(path_new):
                os.remove(path_new)
            os.symlink(path_orig, path_new)
            paths_new.append(path_new)

        # Return output dictionary
        return [{'name': os.path.basename(p), 'path': p} for p in paths_new]

    def __construct_config_file(self, input_files: List[Dict[str, str]]) -> str:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        config_data = self.get_template_data(input_files)

        with CONFIG_DATA.open() as handle_in:
            config_data.update(yaml.load(handle_in.read(), Loader=yaml.SafeLoader))

        return SnakePipelineUtils.generate_config_file(config_data, self._working_dir)

    def get_template_data(self, input_data: [List[Dict[str, str]]]) -> Dict[str, Any]:
        """
        Returns the template data that is common to all pipelines
        :return: Template data
        """
        template_data = {
            'pipeline': {
                'name': self._name,
                'version': f"{self._version}",
            },
            'ubams'      : [("").join(file['name'].split(".")[:-2]) for file in input_data],
            'sample': self._args.sample,
            'ubam_prefix': self._args.ubam_prefix,
            'working_dir': str(self._working_dir),
        }
        return template_data

    def _run_snakemake_main(self, config_file: str) -> None:
        """
        Runs the main snakefile for the pipeline.
        :param config_file: Configuration file
        :return: None
        """
        #if self._pipeline.keep_config:
        #    self._pipeline.log_config_file(config_file)

        try:
            SnakePipelineUtils.run_snakemake(self._snakefile, config_file, [], self._working_dir, self._args.threads)
            log_file = self._working_dir / 'camel.log'
            if log_file.exists():
                shutil.copyfile(str(log_file), str(Path(self._args.output_dir) / 'camel.log'))
            logging.info("Pipeline finished successfully")
        except SnakemakeExecutionError as err:
            if self._pipeline.keep_error_log:
                self._pipeline.log_error_to_file(err)
            raise err


    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser(description='Run the GATK best practices pipeline.')

        # input
        parser.add_argument('--working-dir', dest = "working_dir", default=os.path.abspath('.'), type=str, help='Working directory')
        parser.add_argument('--ubam-prefix', dest = "ubam_prefix", type=str, help='Path to uBAM files')
        parser.add_argument('--sample', dest = "sample", type=str, help='Sample name')

        # logs
        parser.add_argument('--log', action='store_true', help="If this flag is set, config file and error logs are kept")

        # other
        parser.add_argument('--threads', dest = "threads", type=str)


        return parser.parse_args(args)

    @property
    def name(self) -> str:
        """
        Returns the pipeline name.
        :return: Name
        """
        return self._name

    @property
    def version(self) -> str:
        """
        Returns the pipeline version.
        :return: Version
        """
        return self._version

if __name__ == '__main__':
    Camel.get_instance()
    main = MainBroadWGSPipeline()
    main.run()