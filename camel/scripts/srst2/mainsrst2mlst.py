#!/usr/bin/env python
import argparse
import logging
from pathlib import Path
from typing import Optional, Dict, List, Sequence

import shutil

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
from camel.app.tools.srst2.srst2mlst import Srst2Mlst


class MainSrst2Mlst(object):
    """
    Main script to run the base SRST2 sequence typing tool from galaxy.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the wrapper.
        :param args: (Optional) arguments, if not specified arguments are parsed from command line
        """
        self._args = MainSrst2Mlst._parse_arguments(args)
        self._dir_working = Path(self._args.working_dir)
        self._camel = Camel()

    def run(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        srst2 = Srst2Mlst(self._camel)
        srst2.add_input_files(self.__get_reads_input())
        srst2.add_input_files(self.__get_database_input())
        self._add_advanced_options(srst2)
        srst2.run(self._args.working_dir)
        shutil.copyfile(srst2.tool_outputs['TSV'][0].path, self._args.output_tsv)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--fastq-pe', nargs=2, type=str)
        argument_parser.add_argument('--fastq-se', type=str)
        argument_parser.add_argument('--profiles-tsv', type=str,
                                     help="tabular file containing the sequence type definitions")
        argument_parser.add_argument('--locus-fasta', type=str, help="FASTA file containing the sequence for all loci.")
        argument_parser.add_argument('--delimiter', type=str, default='_')
        argument_parser.add_argument('--output-tsv', type=str)
        argument_parser.add_argument('--max-mismatch', type=int)
        argument_parser.add_argument('--min-depth', type=int)
        argument_parser.add_argument('--min-edge-depth', type=int)
        argument_parser.add_argument('--max-unaligned-overlap', type=int)
        argument_parser.add_argument('--working-dir', type=str, default=str(Path('.').absolute()))
        return argument_parser.parse_args(args)

    def __get_reads_input(self) -> Dict[str, List[ToolIOFile]]:
        """
        Returns the reads input.
        Symbolic links are created so SRST2 can recognize forward (_1) & reverse reads (_2) for the paired end .
        :return: Input dictionary
        """
        if self._args.fastq_pe:
            symlink_paths = [self._dir_working / f'reads_{ori}.fastq' for ori in (1, 2)]
            for orig_path, symlink_path in zip(self._args.fastq_pe, symlink_paths):
                symlink_path.symlink_to(orig_path)
            return {'FASTQ_PE': [ToolIOFile(p) for p in symlink_paths]}
        else:
            symlink_path = self._dir_working / 'reads.fastq'
            symlink_path.symlink_to(self._args.fastq_se)
            return {'FASTQ_SE': [ToolIOFile(symlink_path)]}

    def __get_database_input(self) -> Dict[str, List[ToolIOFile]]:
        """
        Returns a dictionary with the input files for the database.
        :return: Dictionary containing the sequences file & the profiles file
        """
        db_path = self._dir_working / 'sequences.fasta'
        db_path.symlink_to(self._args.locus_fasta)
        samtools_faidx = SamtoolsFastaIndex(self._camel)
        samtools_faidx.add_input_files({'FASTA': [ToolIOFile(db_path)]})
        samtools_faidx.run(self._args.working_dir)
        input_dict = {'FASTA': [ToolIOFile(str(db_path))]}
        if self._args.profiles_tsv is not None:
            input_dict['TSV'] = [ToolIOFile(self._args.profiles_tsv)]
        return input_dict

    def _add_advanced_options(self, srst2: Srst2Mlst) -> None:
        """
        Adds the advanced options to the tool.
        :param srst2: SRST2 instance
        :return: None
        """
        if self._args.max_mismatch:
            srst2.update_parameters(max_mismatch=self._args.max_mismatch)
        if self._args.min_depth:
            srst2.update_parameters(min_depth=self._args.min_depth)
        if self._args.min_edge_depth:
            srst2.update_parameters(min_edge_depth=self._args.min_edge_depth)
        if self._args.max_unaligned_overlap:
            srst2.update_parameters(max_unaligned_overlap=self._args.max_unaligned_overlap)


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main_srst2 = MainSrst2Mlst()
    main_srst2.run()
