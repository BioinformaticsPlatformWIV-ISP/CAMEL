#!/usr/bin/env python
import argparse
import logging
import shutil
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from Bio.Seq import Seq

from camel.app.camel import Camel
from camel.app.error.snakemakeexecutionerror import SnakemakeExecutionError
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.screened2 import SNAKEFILE_MAIN


class MainScreened2(object):
    """
    Main class to run the Primer-Probe Check Tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        :return: None
        """
        self._args = MainScreened2._parse_arguments(args)
        self._working_dir = self._args.working_dir
        self._snakefile = SNAKEFILE_MAIN

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return 'Primer-Probe Check Tool'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments (optional)
        :return: Arguments
        """
        parser = argparse.ArgumentParser()

        # Input

        parser.add_argument('--fasta', '-f', type=Path,
                            help='Input Fasta file', required=True)
        parser.add_argument('--primers', '-p', type=Path,
                            help='Tab-delimited config file for primer and probe sequences. '
                                 'Each line should contain in the first column the name of the primer or probe '
                                 'and in the second column the sequence. '
                                 'It is important that in the name of forward primers FW is mentioned, '
                                 'in the name of reverse primers RV and in the name of the probes PR. '
                                 'The names of the primers and probes should not contain spaces.',
                            required=True)
        parser.add_argument('--perc-mismatch', '-pm', type=float,
                            help='Floating point number between 0 and 1 that gives the maximum percentage '
                                 'of allowed mismatches. The default is 0, which is a perfect match',
                            default=0)
        parser.add_argument('--end-mismatch', '-em', type=int,
                            help='Integer between 0 and the total fragment length of the smallest method fragment '
                                 'that should be considered to check for mismatches at the 3 end for the forward '
                                 'and reverse primer. This does NOT apply to the probe.')
        parser.add_argument('--threads', '-t', type=int, help='Number of threads', default=1)
        parser.add_argument('--revcompl-all', '-rca', action="store_true",
                            help='Check the reverse complements of all primers')
        parser.add_argument('--revcompl-list', '-rcl', type=str,
                            help='Check the reverse complement of specific primers (comma-separated)')
        parser.add_argument('--split-fasta-file', '-s', type=float,
                            help='Split fasta file for speed (default: 20)',
                            default=20)

        # Output
        parser.add_argument('--output-dir', '-o', type=Path, help='Output directory', required=True)
        parser.add_argument('--working-dir', '-w', type=Path, help='Working directory', required=True)

        return parser.parse_args(args)

    def _symlink_input(self) -> Dict[str, Any]:
        """
        Symlinks the input files.
        :return: Dictionary with FASTA input
        """
        if not self._working_dir.exists():
            self._working_dir.mkdir(parents=True)
        path_link = self._working_dir / self._args.fasta.name
        if not path_link.is_symlink():
            path_link.symlink_to(self._args.fasta)
        else:
            logging.warning(f'{path_link} exists, not creating link')
        return {'fasta': str(path_link)}

    def _check_perc_mismatch(self) -> float:
        """
        Check if the perc_mismatch argument is a float
        :return: perc_mismatch argument (float)
        """
        value = float(self._args.perc_mismatch)
        if not (0 <= value <= 1):
            raise argparse.ArgumentTypeError(f"{value} is not in the range [0, 1]")
        return value

    def _check_end_mismatch(self) -> int:
        """
        Process primers and check if reverse complements are desired
        :return: end_mismatch argument (int)
        """
        with open(self._args.primers) as handle:
            primer_seq_name_dict = {}
            for line in handle.readlines():
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    primer_name, primer_seq = parts
                    primer_seq_name_dict[primer_name] = [primer_seq]

            # Find the shortest string within each inner list
            shortest_strings = [min(primer_seq_list, key=len) for primer_seq_list in primer_seq_name_dict.values()]

            # Find the overall shortest string
            overall_shortest_string_len = len(min(shortest_strings, key=len))
            if overall_shortest_string_len < int(self._args.end_mismatch):
                self._args.end_mismatch = overall_shortest_string_len
                logging.warning(f'end_mismatch parameter is too big, length of shortest primer is taken')
        return int(self._args.end_mismatch)

    def _primer_dictionary(self) -> Dict[str, Any]:
        """
        Process primers and check if reverse complements are desired
        :return: Dictionary with primers input
        """
        with open(self._args.primers) as primer_file:
            primer_seq_name_dict = {
                primer_name: [primer_seq]
                for line in primer_file
                if (parts := line.strip().split('\t')) and len(parts) == 2
                for primer_name, primer_seq in [parts]
            }

            if self._args.revcompl_all:
                primer_seq_name_dict.update({
                    primer_name + "§RevCompl": [str(Seq(primer_seq[0]).reverse_complement())]
                    for primer_name, primer_seq in primer_seq_name_dict.items()
                })
            elif self._args.revcompl_list:
                rev_complement_path_split = self._args.revcompl_list.split(",")
                primer_seq_name_dict.update({
                    primer_name + "§RevCompl": [str(Seq(primer_seq[0]).reverse_complement())]
                    for primer_name, primer_seq in primer_seq_name_dict.items()
                    if primer_name in rev_complement_path_split
                })
            else:
                logging.warning('Not checking the reverse complements')
        return primer_seq_name_dict

    def _run_snakemake_main(self, config_file: str) -> None:
        """
        Runs the main snakefile for the pipeline.
        :param config_file: Configuration file
        :return: None
        """
        log_file = self._args.working_dir / 'camel.log'
        try:
            # Run snakemake
            SnakePipelineUtils.run_snakemake(
                self._snakefile, config_file, [], self._args.working_dir, self._args.threads)
            logging.info("Pipeline finished successfully")
        except SnakemakeExecutionError:
            raise RuntimeError(f"Error executing Snakemake. Check log for more information: {log_file}")

        # Copy log file to output directory if that directory is given
        if log_file.exists() and 'output_dir' in self._args:
            shutil.copyfile(str(log_file), str(Path(self._args.output_dir) / 'camel.log'))

    def __construct_config_file(self, input_files: Dict[str, str]) -> str:
        """
        Constructs the configuration file.
        :param input_files: Dictionary with fasta file
        :return: Configuration file
        """
        config_data = {
            'input': input_files,
            'primers': self._primer_dictionary(),
            'perc_mismatch': float(self._check_perc_mismatch()),
            'end_mismatch': int(self._check_end_mismatch()),
            'output_dir': str(self._args.output_dir),
            'split_fasta_file': int(self._args.split_fasta_file),
            'working_dir': str(self._args.working_dir),
            'threads': int(self._args.threads),
            'parts_size': int(len([1 for line in open(Path(self._args.fasta)) if
                                   line.startswith(">")]) / int(self._args.split_fasta_file))
        }
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainScreened2()
    main.run()
