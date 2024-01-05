#!/usr/bin/env python
import argparse
import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Sequence, Any

import yaml
from camel.app.camel import Camel
from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.ncbihumanreadscrubber import CONFIG_DATA
from camel.scripts.ncbihumanreadscrubber import SNAKEFILE_MAIN


class MainNcbiHumanReadScrubber(ReportPipeline):
    """
    Main class to run the Ncbi human read scrubber tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('NCBI human read scrubber stand alone', '0.2', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return 'NCBI human read scrubber'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)

    def __construct_config_file(self, input_files: List[Dict[str, str]]) -> str:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        if self._args.fastq_pe is not None:
            key = 'fastq_pe'
        elif self._args.fastq_se is not None:
            key = 'fastq_se'
        else:
            key = 'fasta'
        config_data = self.get_template_data(key, input_files)
        config_data['read_type'] = self._args.read_type
        # Add existing config data
        with open(CONFIG_DATA) as handle_in:
            config_data.update(yaml.safe_load(handle_in.read()))
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    def _get_fastq_input_links(self) -> List[List[Tuple[Path, str]]]:
        """
        Returns the links to the input FASTQ files.
        :return: Links
        """
        links = []
        if self._args.fastq_se is not None:
            gzipped = FileSystemHelper.is_gzipped(self._args.fastq_se)
            links.append([self._args.fastq_se, f"{self.sample_name}.fastq{'.gz' if gzipped else ''}"])
        else:
            for read_nb, path in enumerate(self._args.fastq_pe, start=1):
                gzipped = FileSystemHelper.is_gzipped(path)
                links.append([path, f"{self.sample_name}_{read_nb}.fastq{'.gz' if gzipped else ''}"])
        return links

    @property
    def sample_name(self) -> str:
        """
        Returns the sample name.
        :return: Sample name
        """
        if self._args.fastq_pe is not None:
            return super().sample_name
        elif self._args.fastq_se is not None:
            name = self._args.fastq_se_name if (self._args.fastq_se_name is not None) else self._args.fastq_se
            return FastqUtils.get_sample_name(name, FastqUtils.PATTERN_FQ_SE)
        else:
            # three options for fasta: fasta name, sample name or nothing
            pattern_fasta = r'(.+?)(_S\d+)?(_L\d{3})?(_\d+)?.(fasta|fa)(.gz)?'
            if self._args.sample_name is not None:
                return FileSystemHelper.make_valid(self._args.sample_name)
            elif self._args.fasta_name is not None:
                return FastqUtils.get_sample_name(self._args.fasta_name, pattern_fasta)
            else:
                return FastqUtils.get_sample_name(self._args.fasta, pattern_fasta)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(argument_parser)
        argument_parser.add_argument('--fastq-se', type=Path, help="Input SE FASTQ file")
        argument_parser.add_argument('--fastq-se-name', help="Input SE FASTQ file name")
        argument_parser.add_argument('--fasta', help="Input FASTA file", type=Path)
        argument_parser.add_argument('--fasta-name', help="Name of the input FASTA file", type=str)
        return argument_parser.parse_args(args)

    def _symlink_input(self) -> List[Dict[str, Any]]:
        """
        Symlinks the input files.
        :return: List of FASTQ/FASTA input dictionaries
        """
        # Determine link names if FASTQ input
        if self._args.fastq_pe or self._args.fastq_se:
            links = self._get_fastq_input_links()
        else:
            links = [(self._args.fasta, f"{self.sample_name}.fasta")]

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
            # ToolIO fasta file if fasta
            if self._args.fasta:
                SnakemakeUtils.dump_object(
                    [ToolIOFile(path_new)], self._args.working_dir / 'input' / f"{self.sample_name}.io")
        # Return output dictionary
        return [{'name': p.name, 'path': str(p)} for p in paths_new]


if __name__ == '__main__':
    Camel.get_instance()
    main = MainNcbiHumanReadScrubber()
    main.run()
