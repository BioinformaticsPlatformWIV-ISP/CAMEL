#!/usr/bin/env python
import argparse
import logging
from pathlib import Path
from typing import Optional, Dict, List, Sequence

import shutil

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.srst2.srst2gene import Srst2Gene


class MainSrst2Gene(object):
    """
    Main script to run the base SRST2 gene detection tool from galaxy.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the wrapper.
        :param args: (Optional) arguments, if not specified arguments are parsed from command line
        """
        self._args = MainSrst2Gene._parse_arguments(args)
        self._camel = Camel()

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--fastq-pe', nargs=2, type=str)
        argument_parser.add_argument('--fastq-se', type=Path)
        argument_parser.add_argument('--gene-fasta', type=Path)
        argument_parser.add_argument('--gene-fasta-name', type=str)
        argument_parser.add_argument('--output-tsv', type=Path)
        argument_parser.add_argument('--min-coverage', type=int, default=90)
        argument_parser.add_argument('--max-divergence', type=int, default=10)
        argument_parser.add_argument('--min-depth', type=int)
        argument_parser.add_argument('--min-edge-depth', type=int)
        argument_parser.add_argument('--max-unaligned-overlap', type=int)
        argument_parser.add_argument('--working-dir', type=Path, default=Path.cwd())
        return argument_parser.parse_args(args)

    def __get_reads_input(self) -> Dict[str, List[ToolIOFile]]:
        """
        Returns the reads input.
        Symbolic links are created so SRST2 can recognize forward (_1) & reverse reads (_2) for the paired end .
        :return: Input dictionary
        """
        if self._args.fastq_pe:
            symlink_paths = [self._args.working_dir / f'reads_{ori}.fastq' for ori in (1, 2)]
            for orig_path, symlink_path in zip(self._args.fastq_pe, symlink_paths):
                symlink_path.symlink_to(orig_path)
            return {'FASTQ_PE': [ToolIOFile(p) for p in symlink_paths]}
        else:
            symlink_path = self._args.working_dir / 'reads.fastq'
            symlink_path.symlink_to(self._args.fastq_se)
            return {'FASTQ_SE': [ToolIOFile(symlink_path)]}

    def __add_advanced_options(self, srst2: Srst2Gene) -> None:
        """
        Adds the advanced options.
        :param srst2: SRST2 instance
        :return: None
        """
        if self._args.min_coverage:
            srst2.update_parameters(min_coverage=self._args.min_coverage)
        if self._args.max_divergence:
            srst2.update_parameters(max_divergence=self._args.max_divergence)
        if self._args.min_depth:
            srst2.update_parameters(min_depth=self._args.min_depth)
        if self._args.min_edge_depth:
            srst2.update_parameters(min_edge_depth=self._args.min_edge_depth)
        if self._args.max_unaligned_overlap:
            srst2.update_parameters(max_unaligned_overlap=self._args.max_unaligned_overlap)

    def run(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        srst2 = Srst2Gene(self._camel)

        # Add input files
        srst2.add_input_files(self.__get_reads_input())
        db_name = Path(self._args.gene_fasta_name if self._args.gene_fasta_name is not None else self._args.gene_fasta)\
            .stem
        db_path = self._args.working_dir / db_name
        db_path.symlink_to(self._args.gene_fasta)
        srst2.add_input_files({'FASTA': [ToolIOFile(db_path)]})

        # Add options
        self.__add_advanced_options(srst2)

        # Run
        srst2.run(self._args.working_dir)

        # Collect output
        if 'TSV' in srst2.tool_outputs:
            shutil.copyfile(srst2.tool_outputs['TSV'][0].path, self._args.output_tsv)
        else:
            logging.info("No genes detected")


if __name__ == '__main__':
    Camel.get_instance()
    srst2_Wrapper = MainSrst2Gene()
    srst2_Wrapper.run()
