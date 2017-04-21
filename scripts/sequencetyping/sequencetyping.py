import argparse
import logging
import os
import shutil
import traceback

from app.camel import Camel
from app.components.files.fastqutils import FastqUtils
from app.io.tooliodirectory import ToolIODirectory
from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue
from app.pipeline.pipeline import Pipeline
from config import DB_CONFIG, LOGGING_CONFIG
from resources import YAML_TYPING_FAST, YAML_TRIMMING_ASSEMBLY_SPADES
from scripts.sequencetyping.htmlreporter import HtmlReporter


class MainMlst(object):
    """
    Class to run sequence typing workflow.
    """

    def __init__(self):
        """
        Initializes the workflow.
        """
        self._args = MainMlst._parse_arguments()
        self._camel = Camel(DB_CONFIG, LOGGING_CONFIG)
        self._report = HtmlReporter(self._args.output_dir, self._args.output_html)
        self._pipeline_assembly = None

    @staticmethod
    def _parse_arguments():
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--scheme-dir', required=True)
        argument_parser.add_argument('--fasta')
        argument_parser.add_argument('--fasta-name')
        argument_parser.add_argument('--fastq', nargs=2)
        argument_parser.add_argument('--fastq-names', nargs=2)
        argument_parser.add_argument('--output-dir', required=True)
        argument_parser.add_argument('--output-html', required=True)
        argument_parser.add_argument('--kmers')
        return argument_parser.parse_args()

    def run(self):
        """
        Runs the workflow.
        :return: None
        """
        self._report.initialize(self._args.fasta_name, self._args.fastq_names)
        try:
            if self._args.fastq:
                assembly = self.__run_assembly_pipeline()
            else:
                assembly = ToolIOFile(self._args.fasta)
            self.__run_sequence_typing_pipeline(assembly)
        except ValueError as err:
            logging.error(traceback.format_exc())
            raise err
        finally:
            shutil.copy('camel.log', os.path.join(self._args.output_dir, 'log.txt'))

    def __run_assembly_pipeline(self):
        """
        Runs the de-novo assembly pipeline.
        :return: Assembly output file object
        """
        pipeline_assembly = Pipeline([YAML_TRIMMING_ASSEMBLY_SPADES], self._camel)
        pipeline_assembly.set_initial_input({'FASTQ_PE': [ToolIOFile(fq_file) for fq_file in self._args.fastq]})
        if self._args.kmers:
            pipeline_assembly.add_job_options({'Assembly': {'kmers': self._args.kmers}})
        pipeline_assembly.run('output')
        self._report.add_trimming_section(pipeline_assembly.get_step('Read_trimming').informs)
        self._report.add_fastqc_section(pipeline_assembly.get_step('FastQC').outputs['HTML'])
        self._report.add_assembly_section(
            self.__get_sample_name(), pipeline_assembly.get_step('Assembly').outputs['FASTA_Contig'][0])
        return pipeline_assembly.get_step('Assembly').outputs['FASTA_Contig'][0]

    def __get_sample_name(self):
        """
        Returns the sample name.
        :return: None
        """
        if self._args.fasta_name is not None:
            return os.path.splitext(self._args.fasta_name)[0]
        else:
            try:
                return FastqUtils.get_sample_name(self._args.fastq_names[0])
            except ValueError:
                logging.warning('Cannot determine sample name from FASTQ names, using default name')
        return 'na'

    def __run_sequence_typing_pipeline(self, assembly):
        """
        Runs the sequence typing pipeline.
        :param assembly: Assembly
        :return: None
        """
        pipeline = Pipeline([YAML_TYPING_FAST], self._camel, True)
        pipeline.set_initial_input({
            'DIR': [ToolIODirectory(self._args.scheme_dir)],
            'FASTA': [assembly],
            'DIR_HTML': [ToolIODirectory(self._args.output_dir)],
            'HTML': [ToolIOFile(self._args.output_html)],
            'SAMPLE_NAME': [ToolIOValue(self.__get_sample_name())]
        })
        pipeline.run('output')

if __name__ == '__main__':
    sequence_typing = MainMlst()
    sequence_typing.run()
