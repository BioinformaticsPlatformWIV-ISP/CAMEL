#!/usr/bin/env python
import argparse
import glob
import logging
import os
import shutil
import yaml
from pathlib import Path
from typing import Optional, Sequence

from camel.app.camel import Camel
from camel.app.error.snakemakeexecutionerror import SnakemakeExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.pipeline import Pipeline
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.broadwgs import SNAKEFILE_MAIN, CONFIG_DATA, REFERENCES, INTERVALS, TOOL_DATA


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
        self._version = '1.0'
        self._snakefile = SNAKEFILE_MAIN
        self._args = MainBroadWGSPipeline._parse_arguments(args)
        self._working_dir = Path(self._args.working_dir)
        self._pipeline = Pipeline(self._name, Camel.get_instance(), self._args.log, self._args.log)

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        self._symlink_input()

        config_file = self.__construct_config_file()

        try:
            resources = dict([arg.split(",") for arg in self._args.resources])
            SnakePipelineUtils.run_snakemake(self._snakefile, config_file, [], self._working_dir, self._args.threads, resources)
            log_file = self._working_dir / 'camel.log'
            if log_file.exists():
                shutil.copyfile(str(log_file), str(Path(self._args.working_dir) / 'output' / 'camel.log'))
            logging.info("Pipeline finished successfully")
        except SnakemakeExecutionError as err:
            if self._pipeline.keep_error_log:
                self._pipeline.log_error_to_file(err)
            raise err

    def _symlink_input(self) -> None:
        """
        Symlinks the input files.
        :return: None
        """
        base_path = self._args.input

        # Create directory
        dir_links = self._working_dir / 'input'
        if not dir_links.exists():
            dir_links.mkdir(parents=True, exist_ok=True)

        # Determine link names + generate ToolIOFiles
        links = []
        for file in base_path:
            #add extension
            r1 = file + "R1.fastq.gz"
            r2 = file + "R2.fastq.gz"
            #append to links for symlinking
            links.append(r1)
            links.append(r2)
            #ToolIOFile
            io_files = [ToolIOFile(r1), ToolIOFile(r2)]
            SnakemakeUtils.dump_object(io_files, os.path.join(dir_links, os.path.basename(file) + ".fastq.gz.io"))

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

    def __construct_config_file(self) -> str:
        """
        Constructs the configuration file
        :return: Configuration file
        """
        # pipeline information and cmd line arguments

        config_data = {
            'pipeline': {
                'name': self._name,
                'version': f"{self._version}",
            },
            'params_smk': {
                'threads': self._args.threads,
                'threads_bwa': round(self._args.threads * 0.9 / len(self._args.input)),
                'resources': self._args.resources,
            },
            'sample': self._args.sample,
            'input_basenames': [os.path.basename(f) for f in self._args.input],
            'working_dir': str(self._working_dir),
            'final_output_dir': str(self._working_dir / 'output')
        }

        # add data from config yml
        with CONFIG_DATA.open() as handle_in:
            config_data.update(yaml.load(handle_in.read(), Loader=yaml.SafeLoader))

        # add data from user specified references yml if defined
        if self._args.references is not None:
            reference_file = self._args.references
        else:
            reference_file = REFERENCES
        with reference_file.open() as handle_in:
            config_data.update(yaml.load(handle_in.read(), Loader=yaml.SafeLoader))

        # add intervals - based on interval files
        intervals_list = []
        if self._args.intervals is not None:
            intervals_location = self._args.intervals
        else:
            intervals_location = INTERVALS
        interval_files = glob.glob(os.path.join(intervals_location, "interval_*.intervals"))
        config_data.update({"intervals": list(range(len(interval_files))),
                            "intervals_location": str(intervals_location)})

        # add tool parameters from tool_data yml
        if self._args.tool_data is not None:
            tool_data = self._args.tool_data
        else:
            tool_data = TOOL_DATA
        with tool_data.open() as handle_in:
            config_data.update(yaml.load(handle_in.read(), Loader=yaml.SafeLoader))

        # Update read length
        config_data["rule_params"]["qc"]["picard_wgs_metrics"]["read_length"] = self._args.read_length
        config_data["rule_params"]["qc"]["picard_raw_wgs_metrics"]["read_length"] = self._args.read_length

        return SnakePipelineUtils.generate_config_file(config_data, self._working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser(description='Run the GATK best practices pipeline.')

        # input data
        parser.add_argument('-i', dest = "input", nargs="+", help="FastQ input files")
        parser.add_argument('--working-dir', dest = "working_dir", default = str(Path('.').absolute()), type=str, help='Working directory')
        parser.add_argument('--sample', type = str, help = 'Sample name')
        parser.add_argument('--read-length', dest = "read_length", type = int, default = 250, help = 'Read length. Default: 250')
        parser.add_argument('--references', type = str, help = "Path to references yml")
        parser.add_argument('--intervals', type = str, help = 'Path to interval files directory')
        parser.add_argument('--tool-data', dest = "tool_data", type = str, help = "Path to tool data yml")

        # logs
        parser.add_argument('--log', action='store_true', help="If this flag is set, config file and error logs are kept")

        # Snakemake parameters
        parser.add_argument('--threads', type=int, help="Snakemake parameter: number of cores")
        parser.add_argument('--resources', type=str, nargs="+", help="Snakemake parameter: resources. Key-value pair separated by a comma, e.g. mem_mb,1000. Multiple pairs allowed")


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