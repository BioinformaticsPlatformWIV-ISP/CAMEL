import argparse
import datetime
import logging
from dataclasses import dataclass
from typing import List, Optional, Dict

import os
import shutil

from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.galaxy.galaxyutils import GalaxyUtils
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.workflows.assemblywrapper import AssemblyWrapper
from camel.app.components.workflows.readtrimmingwrapper import ReadTrimmingWrapper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources import CSS_STYLE
from camel.resources.javascript import JQUERY_SRC


@dataclass
class ReadInput:
    pe: List[ToolIOFile]
    se_fwd: List[ToolIOFile]
    se_rev: List[ToolIOFile]


class MainScriptHelper(object):
    """
    This class is used as helper class for common operation for main scripts.
    For example the read trimming, de-novo assembly that is shared by tools such as 'MLST tool', 'Gene detection tool',
    'ResFinder local' can be done by this class.
    """

    def __init__(self, working_dir: str, sample_name: str):
        """
        Initializes the helper class.
        :param working_dir: Working directory
        """
        self._working_dir = working_dir
        self._sample_name = sample_name
        self._log_files = {}
        self._informs = []

    def trim_reads(self, fastq_reads_raw: List[str], report: Optional[HtmlReport] = None, threads: int = 8,
                   export_fastq: bool = False) -> ReadInput:
        """
        Runs the read trimming workflow.
        :param fastq_reads_raw: Raw FASTQ PE reads
        :param report: If set, the output is added to the given report
        :param threads: Number of threads
        :param export_fastq: If True, FASTQ files are included in the report
        :return: Tuple of trimmed PE reads, trimmed forward reads, trimmed reverse reads
        """
        logging.info("Trimming FASTQ input reads")
        trimming = ReadTrimmingWrapper(os.path.join(self._working_dir, 'trimming'))
        trimming.run_workflow(fastq_reads_raw, threads, export_fastq)
        if report is not None:
            report.add_html_object(trimming.output.report_section)
            trimming.output.report_section.copy_files(report.output_dir)
            report.save()
        if trimming.output.log_file is not None:
            self._log_files['trimming'] = trimming.output.log_file
        self._informs.append(trimming.output.informs_trimmomatic)
        return ReadInput(trimming.output.trimmed_reads_pe, trimming.output.trimmed_reads_se_fwd,
                         trimming.output.trimmed_reads_se_rev)

    def assemble_fastq_reads(self, assembly_input: ReadInput, report: Optional[HtmlReport] = None,
                             kmers: Optional[str] = None, min_contig_length: int = 0, cov_cutoff: str = 'off',
                             threads: int = 8) -> ToolIOFile:
        """
        Assembles FASTQ reads using SPAdes
        :param assembly_input: Assembly input
        :param report: If set, the output is added to the given report
        :param kmers: Comma separated list of Kmers to use for the assembly
        :param min_contig_length: Minimal contig length
        :param cov_cutoff: Contig coverage cutoff
        :param threads: Number of threads to use
        :return: ToolIOFile FASTA object with the assembled contigs
        """
        logging.info("Starting de-novo assembly")
        assembly = AssemblyWrapper(os.path.join(self._working_dir, 'assembly'))
        assembly.run_workflow(
            self._sample_name, assembly_input.pe, assembly_input.se_fwd, assembly_input.se_rev, kmers, cov_cutoff,
            min_contig_length, threads)
        if report is not None:
            report.add_html_object(assembly.output.report_section)
            assembly.output.report_section.copy_files(report.output_dir)
            report.save()
        if assembly.output.log_file is not None:
            self._log_files['assembly'] = assembly.output.log_file
        self._informs.extend(assembly.output.informs)
        return assembly.output.fasta_contigs

    @staticmethod
    def determine_sample_name(args: argparse.Namespace) -> str:
        """
        Determines the sample names based on the given command line arguments.
        :return: Sample name
        """
        if ('sample_name' in args) and (args.sample_name is not None):
            return args.sample_name
        elif ('fasta_name' in args) and (args.fasta_name is not None):
            return os.path.splitext(args.fasta_name)[0]
        elif ('fasta' in args) and (args.fasta is not None):
            return os.path.splitext(os.path.basename(args.fasta))[0]
        elif (args.fastq_pe_names is not None) or (args.fastq_pe is not None):
            names = args.fastq_pe_names if args.fastq_pe_names else [os.path.basename(x) for x in args.fastq_pe]
            try:
                # See if it matches a standard FASTQ format
                return FastqUtils.get_sample_name(names[0])
            except ValueError:
                # Check if it matches a Galaxy format
                return GalaxyUtils.determine_sample_name_from_fq(names, 'NA')
        logging.warning("Cannot determine sample name from given arguments")
        return 'NA'

    @staticmethod
    def determine_input_files(args: argparse.Namespace) -> str:
        """
        Determines the input files based on the given command line arguments.
        :param args: Command line arguments
        :return: Input files as string
        """
        if (args.fasta is not None) and (args.fasta_name is not None):
            return args.fasta_name
        elif args.fasta is not None:
            return os.path.basename(args.fasta)
        elif (args.fastq_pe is not None) and (args.fastq_pe_names is not None):
            return ', '.join(args.fastq_pe_names)
        elif args.fastq_pe is not None:
            return ', '.join([os.path.basename(f) for f in args.fastq_pe])
        logging.warning("Cannot determine input files from given arguments")
        return 'NA'

    @staticmethod
    def add_common_arguments(argument_parser: argparse.ArgumentParser) -> None:
        """
        Adds the common arguments to the argument parser.
        :param argument_parser: Argument parser
        :return: None
        """
        argument_parser.add_argument('--sample-name', type=str)
        argument_parser.add_argument('--output-dir', required=True, type=str)
        argument_parser.add_argument('--output-html', required=True, type=str)
        argument_parser.add_argument('--working-dir', default=os.path.abspath('.'), type=str)
        argument_parser.add_argument('--threads', default=8, type=int)

    @staticmethod
    def add_input_files_arguments(argument_parser: argparse.ArgumentParser) -> None:
        """
        Adds the arguments for the input files (FASTA / FASTQ PE).
        :param argument_parser: Argument parser
        :return: None
        """
        argument_parser.add_argument('--fasta', help="Input FASTA file", type=str)
        argument_parser.add_argument('--fasta-name', help="Input FASTA file name", type=str)
        argument_parser.add_argument('--fastq-pe', help="Input PE FASTQ files", nargs=2)
        argument_parser.add_argument('--fastq-pe-names', help="Input PE FASTQ file names", nargs=2)
        argument_parser.add_argument('--trim-reads', help="Perform read trimming", action='store_true')

    @staticmethod
    def add_assembly_arguments(argument_parser: argparse.ArgumentParser) -> None:
        """
        Adds the arguments that are used for the assembly.
        :param argument_parser: Argument parser
        :return: None
        """
        argument_parser.add_argument('--assembly-kmers', help="Kmers to use for assembly", type=str)
        argument_parser.add_argument(
            '--assembly-cov-cutoff', help="Minimal k-mer coverage for assembled contigs", type=int)
        argument_parser.add_argument(
            '--assembly-min-contig-length', help="Minimal length for assembled contigs", type=int)

    @staticmethod
    def prepare_galaxy_output(output_dir: str, output_html: str) -> None:
        """
        Prepares the Galaxy output files at the start of the script.
        - The output HTML file is removed, so Snakemake can regenerate it
        - The output directory is created if it does not exist yet.
        :param output_dir: Output directory
        :param output_html: Output report path
        :return: None
        """
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        if os.path.isfile(output_html):
            os.remove(output_html)

    @property
    def logs(self) -> Dict[str, str]:
        """
        Returns the log files (key: name, value: log file path).
        :return: Logs
        """
        return self._log_files

    @property
    def informs(self) -> List[Dict[str, str]]:
        """
        Returns the informs.
        :return: List of informs
        """
        return self._informs

    def symlink_input_files(self, fasta_input: Optional[str] = None, fastq_pe_input: Optional[List[str]] = None) \
            -> Dict[str, List[str]]:
        """
        Symlinks the input files.
        :param fasta_input: FASTA input
        :param fastq_pe_input: FASTQ PE input
        :return: None
        """
        # Determine link locations
        links = []
        if fasta_input is not None:
            links.append(['fasta', fasta_input, f'{FileSystemHelper.make_valid(self._sample_name)}.fasta'])
        if fastq_pe_input is not None:
            for read_nb, fq in enumerate(fastq_pe_input, start=1):
                gzipped = FileSystemHelper.is_gzipped(fq)
                filename = f"{FileSystemHelper.make_valid(self._sample_name)}_{read_nb}.fastq{'.gz' if gzipped else ''}"
                links.append(['fastq_pe', fq, filename])

        # Create directory
        link_dir = os.path.join(self._working_dir, 'input')
        if not os.path.isdir(link_dir):
            os.makedirs(link_dir)

        # Create symlinks
        input_files = {}
        for key, path_orig, link_name in links:
            path_link = os.path.join(link_dir, link_name)
            if os.path.islink(path_link):
                os.remove(path_link)
            logging.debug(f"Creating symbolic link: '{path_orig}' to '{path_link}'")
            os.symlink(path_orig, path_link)
            try:
                input_files[key].append(path_link)
            except KeyError:
                input_files[key] = [path_link]

        if len(input_files) == 0:
            raise ValueError("No input files found.")
        return input_files

    def get_blast_input(self, input_files: Dict[str, List[str]], args: argparse.Namespace,
                        report: Optional[HtmlReport] = None) -> ToolIOFile:
        """
        Returns the input for BLAST based detection methods.
        Takes into accounts:
        - Type of input (FASTA / FASTQ)
        - Read trimming (optional)
        - Kmer setting for de-novo assembly
        :param input_files: Input files in a standardized format.
        :param args: Command line arguments
        :param report: HTML report
        :return: FASTA input for BLAST.
        """
        # Return FASTA file if there is one
        if 'fasta' in input_files:
            return ToolIOFile(input_files['fasta'][0])

        # Trim reads if it is specified in the arguments
        if args.trim_reads:
            assembly_input = self.trim_reads(input_files['fastq_pe'], report, args.threads, args.report_include_fastq)
        else:
            assembly_input = ReadInput([ToolIOFile(l) for l in input_files['fastq_pe']], [], [])

        # Perform de-novo assembly
        if args.assembly_cov_cutoff is None:
            cov_cutoff = 'off'
        elif args.assembly_cov_cutoff == 0:
            cov_cutoff = 'auto'
        else:
            cov_cutoff = str(args.assembly_cov_cutoff)
        return self.assemble_fastq_reads(
            assembly_input, report, args.assembly_kmers, args.assembly_min_contig_length, cov_cutoff, args.threads)

    def get_srst2_input(self, input_files: Dict[str, List[str]], args: argparse.Namespace,
                        report: Optional[HtmlReport] = None) -> List[ToolIOFile]:
        """
        Returns the input for SRST2 (forward + reverse reads).
        :param input_files: Input files in a standardized format.
        :param args: Command line arguments
        :param report: HTML report
        :return: SRST2 input files
        """
        if args.trim_reads is False:
            return [ToolIOFile(p) for p in input_files['fastq_pe']]
        else:
            trimming_out = self.trim_reads(input_files['fastq_pe'], report, args.threads, args.report_include_fastq)
            return trimming_out.pe

    @staticmethod
    def init_report(output_path: str, output_dir: str, title: str, header: str) -> HtmlReport:
        """
        Initializes the HTML report.
        :param output_path: Output path
        :param output_dir: Output directory
        :param title: Report title
        :param header: Report header
        :return: Report
        """
        report = HtmlReport(output_path, output_dir, [JQUERY_SRC])
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        report.initialize(title, CSS_STYLE)
        report.add_pipeline_header(header)
        report.save()
        return report

    def export_log_files(self, output_dir: str) -> None:
        """
        Exports the log files to the output directory.
        :param output_dir: Output directory
        :return: None
        """
        dir_logs = os.path.join(output_dir, 'logs')
        if not os.path.isdir(dir_logs):
            os.makedirs(dir_logs)
        for key, path in self.logs.items():
            shutil.copyfile(path, os.path.join(dir_logs, f'log_{key}.txt'))

    def export_output_and_commands_section(self, report: HtmlReport, section: HtmlReportSection) -> None:
        """
        Adds the output and commands sections to the report.
        Copies the log files to the output folder.
        :param report: Report
        :param section: Section to add
        :return: None
        """
        report.add_html_object(section)
        section.copy_files(report.output_dir)
        self.export_log_files(report.output_dir)
        if len(self._informs) > 0:
            section_commands = SnakePipelineUtils.create_commands_section(self._informs, self._working_dir)
            report.add_html_object(section_commands)
        report.save()

    @staticmethod
    def export_analysis_info_section(
            report: HtmlReport, input_files_str: str, additional_info: Optional[List[List[str]]] = None) -> None:
        """
        Exports the analysis info section.
        :param report: Report
        :param input_files_str: Input files as a string
        :param additional_info: Additional information for the info table
        :return: None
        """
        section = HtmlReportSection('Analysis info')
        data = [
            ['Analysis date:', datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT)],
            ['Input file(s):', input_files_str],
        ]
        if additional_info is not None:
            data.extend(additional_info)
        section.add_table(data, table_attributes=[('class', 'information')])
        report.add_html_object(section)
        report.save()
