import argparse
import logging
import os
import shutil

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
from camel.app.tools.srst2.srst2mlst import Srst2Mlst
from camel.config import DB_CONFIG, LOGGING_CONFIG


class MainSrst2Mlst(object):
    """
    Main script to run SRST2 MLST.
    """

    def __init__(self):
        """
        Initializes the wrapper.
        """
        self._args = MainSrst2Mlst._parse_arguments()
        self._camel = Camel(DB_CONFIG, LOGGING_CONFIG)

    def run(self):
        """
        Runs this tool.
        :return: None
        """
        srst2 = Srst2Mlst(self._camel)
        srst2.add_input_files(self._get_reads_input())
        if self._args.scheme_dir is not None:
            srst2.add_input_files(self._get_mlst_files_from_dir())
        else:
            srst2.add_input_files(self._get_mlst_files_from_args())
        self._add_advanced_options(srst2)
        srst2.run('.')
        shutil.copyfile(srst2.tool_outputs['TSV'][0].path, self._args.output_tsv)

    @staticmethod
    def _parse_arguments():
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--fastq-pe', nargs=2)
        argument_parser.add_argument('--fastq-se')
        argument_parser.add_argument('--scheme-dir')
        argument_parser.add_argument('--profiles')
        argument_parser.add_argument('--sequences')
        argument_parser.add_argument('--delimiter')
        argument_parser.add_argument('--output-tsv')
        argument_parser.add_argument('--max-mismatch')
        argument_parser.add_argument('--min-depth')
        argument_parser.add_argument('--min-edge-depth')
        argument_parser.add_argument('--max-unaligned-overlap')
        return argument_parser.parse_args()

    def _get_reads_input(self):
        """
        Returns the reads input. Symbolic links are created so SRST2 can recognize forward (_1) & reverse reads (_2).
        :return: Input dictionary
        """
        if self._args.fastq_pe:
            os.symlink(self._args.fastq_pe[0], 'reads_1.fq')
            os.symlink(self._args.fastq_pe[1], 'reads_2.fq')
            return {'FASTQ_PE': [ToolIOFile('reads_1.fq'), ToolIOFile('reads_2.fq')]}
        else:
            return {'FASTQ_SE': [ToolIOFile(self._args.fastq_se)]}

    def _get_mlst_files_from_dir(self):
        """
        Returns the MLST files from the scheme directory.
        :return: Dictionary containing the sequences file & the profiles file
        """
        logging.info("Checking directory: {}".format(self._args.scheme_dir))
        sequences = ToolIOFile(os.path.join(self._args.scheme_dir, 'sequences.fasta'))
        if not sequences.is_valid():
            raise IOError("Cannot find sequences in scheme directory")
        profiles = ToolIOFile(os.path.join(self._args.scheme_dir, 'profiles.tsv'))
        if not sequences.is_valid():
            raise IOError("Cannot find profiles in scheme directory")
        return {'FASTA': [sequences], 'TSV': [profiles]}

    def _get_mlst_files_from_args(self):
        """
        Returns the MLST files provided in the command line arguments.
        :return: Dictionary containing the sequences file & the profiles file
        """
        os.symlink(self._args.sequences, 'input_sequences.fasta')
        samtools_faidx = SamtoolsFastaIndex(self._camel)
        samtools_faidx.add_input_files({'FASTA': [ToolIOFile('input_sequences.fasta')]})
        samtools_faidx.run(os.path.abspath('.'))
        return {'TSV': [ToolIOFile(self._args.profiles)], 'FASTA': [ToolIOFile(self._args.sequences)]}

    def _add_advanced_options(self, srst2):
        """
        Adds the advanced options.
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
    main_srst2 = MainSrst2Mlst()
    main_srst2.run()
