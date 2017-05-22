import argparse
import shutil

import logging
import os

from app.camel import Camel
from app.io.tooliofile import ToolIOFile
from app.tools.srst2.srst2gene import Srst2Gene
from config import DB_CONFIG, LOGGING_CONFIG


class MainSrst2Gene(object):
    """

    Class to run SRST2 from galaxy.
    """

    def __init__(self):
        """
        Initializes the wrapper.
        """
        self._args = MainSrst2Gene._parse_arguments()
        self._camel = Camel(DB_CONFIG, LOGGING_CONFIG)

    @staticmethod
    def _parse_arguments():
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--fastq-pe', nargs=2)
        argument_parser.add_argument('--fastq-se')
        argument_parser.add_argument('--gene-local-db')
        argument_parser.add_argument('--gene-fasta')
        argument_parser.add_argument('--gene-fasta-name')
        argument_parser.add_argument('--output-tsv')
        argument_parser.add_argument('--min-coverage')
        argument_parser.add_argument('--max-divergence')
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

    def _add_advanced_options(self, srst2):
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

    def run(self):
        """
        Runs this tool.
        :return: None
        """
        srst2 = Srst2Gene(self._camel)
        srst2.add_input_files(self._get_reads_input())
        if self._args.gene_local_db is not None:
            srst2.add_input_files({'FASTA': [ToolIOFile(self._args.gene_local_db)]})
        else:
            os.symlink(self._args.gene_fasta, os.path.basename(self._args.gene_fasta_name))
            srst2.add_input_files({'FASTA': [ToolIOFile(os.path.basename(self._args.gene_fasta_name))]})
        self._add_advanced_options(srst2)
        srst2.run('.')
        if 'TSV' in srst2.tool_outputs:
            shutil.copyfile(srst2.tool_outputs['TSV'][0].path, self._args.output_tsv)
        else:
            logging.info("No genes detected")


if __name__ == '__main__':
    srst2_Wrapper = MainSrst2Gene()
    srst2_Wrapper.run()
