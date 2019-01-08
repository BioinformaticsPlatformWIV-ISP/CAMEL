import argparse
import datetime
from dataclasses import dataclass
from typing import List, Dict, Optional, Union

import os
from Bio import SeqIO

from camel.app.camel import Camel
from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.components.phylogeny.megautils import MEGAUtils
from camel.app.components.phylogeny.newickutils import NewickUtils
from camel.app.components.workflows.readtrimmingwrapper import ReadTrimmingWrapper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.mega.mltreeconstruction import MLTreeConstruction
from camel.app.tools.mega.modelselection import ModelSelection
from camel.app.tools.snpmatrix.snpmatrixconstructor import SnpMatrixConstructor
from camel.resources import CSS_STYLE
from camel.scripts.snpphylogeny import SNAKEFILE_TRIMMING_ALL
from camel.scripts.snpphylogeny.snakefile.trimming_all import TRIMMING_ALL


class InvalidInputError(ValueError):
    pass


class SnpPhylogenyUtils(object):
    """
    Utility class for the SNP phylogeny pipelines.
    """

    @dataclass(frozen=True)
    class Sample:
        """
        This class is used to represent an input sample.
        """
        name_full: str
        name_valid: str
        reads_raw: List[ToolIOFile]
        reads_names: List[str]

        def __hash__(self) -> int:
            """
            Returns the hash value for this
            :return: Hash value
            """
            return hash(self.name_valid)

        def __eq__(self, other: 'SnpPhylogenyUtils.Sample') -> bool:
            """
            Checks if two sample instances are equal.
            :param other: Other instance
            :return: True if equal, False otherwise
            """
            return self.name_valid == other.name_valid

    @dataclass
    class MappingInput:
        """
        This class contains the input for the read mapping step.
        """
        pe: List[ToolIOFile]
        se_fwd: Optional[ToolIOFile] = None
        se_rev: Optional[ToolIOFile] = None

        def as_dict(self) -> Dict[str, str]:
            """
            Return the mapping input as a dictionary
            :return: Mapping input in YAML format
            """
            d = {'PE': [x.path for x in self.pe]}
            if self.se_fwd is not None:
                d['SE_FWD'] = self.se_fwd.path
            if self.se_rev is not None:
                d['SE_REV'] = self.se_rev.path
            return d

    @staticmethod
    def initialize_report(pipeline_name: str, args: argparse.Namespace) -> HtmlReport:
        """
        Initializes a HTML report.
        :param pipeline_name: Pipeline name (e.g. 'Samtools', 'PHEnix')
        :param args: Pipeline arguments
        :return: Initializes HTML report
        """
        report = HtmlReport(args.output_html, args.output_dir)
        if not os.path.isdir(args.output_dir):
            os.makedirs(args.output_dir)
        report.initialize(f'SNP Phylogeny - {pipeline_name}', CSS_STYLE)
        report.add_pipeline_header(f'SNP Phylogeny ({pipeline_name})')
        section = HtmlReportSection('Analysis info')
        section.add_table([
            ['Analysis date:', datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT)],
            ['Nb. of samples:', len(args.sample)],
            ['Reference:', args.reference_name]], table_attributes=[('class', 'information')]
        )
        report.add_html_object(section)
        report.save()
        return report

    @staticmethod
    def add_common_arguments(argument_parser: argparse.ArgumentParser) -> None:
        """
        Adds the common arguments to the given argument parser.
        :param argument_parser: Argument parser
        :return: None
        """
        argument_parser.add_argument('--output-dir', required=True, type=str)
        argument_parser.add_argument('--output-html', required=True, type=str)
        argument_parser.add_argument('--working-dir', default=os.path.abspath('.'), type=str)
        argument_parser.add_argument('--threads', default=8, type=int)
        argument_parser.add_argument('--reference', required=True, type=str)
        argument_parser.add_argument('--reference-name', type=str)
        argument_parser.add_argument('--sample', nargs=5, action='append', required=True)
        argument_parser.add_argument('--trim-reads', action='store_true')
        argument_parser.add_argument('--missing-data', choices=[
            'complete_deletion', 'use_all_sites', 'partial_deletion'])
        argument_parser.add_argument('--site-cov-cutoff', choices=range(0, 101), type=int)
        argument_parser.add_argument('--branch-swap', choices=[
            'none', 'weak', 'very_weak', 'moderate', 'strong', 'very_strong'])
        argument_parser.add_argument('--bootstraps', type=int, required=True)
        argument_parser.add_argument('--ml-method', choices=['nni', 'spr3', 'spr5'], required=True)

    @staticmethod
    def extract_samples(args: argparse.Namespace) -> List['SnpPhylogenyUtils.Sample']:
        """
        Returns the sample names (sorted alphabetically).
        :param args: Command line arguments
        :return: Samples
        """
        samples = []
        for sample_name, fwd_name, fwd_read, rev_name, rev_read in args.sample:
            if (sample_name is None) or (len(sample_name) == 0):
                sample_name = FastqUtils.get_sample_name(fwd_name)
            samples.append(SnpPhylogenyUtils.Sample(
                sample_name, FileSystemHelper.make_valid(sample_name),
                [ToolIOFile(fwd_read), ToolIOFile(rev_read)],
                [fwd_name, rev_name]
            ))

        # Check for duplicate names
        if len(set(s.name_valid for s in samples)) != len(args.sample):
            sample_names = [s.name_valid for s in samples]
            duplicate_sample_names = [x for x in set(sample_names) if sample_names.count(x) > 1]
            raise InvalidInputError("Duplicate sample names are not allowed. Conflicting sample(s): {}".format(
                ', '.join(["'{}' ({} times)".format(n, sample_names.count(n)) for n in duplicate_sample_names])))

        # Check for empty sample names
        empty_samples = [s for s in samples if len(s.name_valid) == 0]
        if len(empty_samples) > 0:
            raise InvalidInputError("Empty sample name for sample(s): {}".format(', '.join([
                str(s.reads_names) for s in empty_samples])))

        # Sort samples by name
        return sorted(samples, key=lambda s: s.name_valid)

    @staticmethod
    def trim_all_reads(samples: List[Sample], working_dir: str, threads: int = 8) -> Dict[
            'SnpPhylogenyUtils.Sample', ReadTrimmingWrapper.ReadTrimmingOutput]:
        """
        Trims all the reads in parallel using Snakemake.
        :param samples: List of samples
        :param working_dir: Working directory
        :param threads: Number of threads
        :return: None
        """
        config_data = {
            'working_dir': working_dir,
            'samples': {s.name_valid: [f.path for f in s.reads_raw] for s in samples}}
        output_file = os.path.join(working_dir, TRIMMING_ALL)
        SnakePipelineUtils.run_snakemake(
            SNAKEFILE_TRIMMING_ALL, config_data, [output_file], working_dir, threads)
        trimming_out_by_sample = SnakemakeUtils.load_object(output_file)
        return {s: trimming_out_by_sample[s.name_valid] for s in samples}

    @staticmethod
    def add_trimming_section_empty(report: HtmlReport) -> None:
        """
        Adds an empty trimming section to the report.
        :param report: Report
        :return: None
        """
        section = HtmlReportSection('Read trimming')
        section.add_paragraph('Read trimming disabled.')
        report.add_html_object(section)
        report.save()

    @staticmethod
    def add_trimming_section(report: HtmlReport, trimming_output_by_sample: Dict[
            'SnpPhylogenyUtils.Sample', ReadTrimmingWrapper.ReadTrimmingOutput]) -> None:
        """
        Adds the trimming section to the report.
        :param report: HTML report
        :param trimming_output_by_sample: Trimming output dictionary
        :return: None
        """
        section = HtmlReportSection('Read trimming')

        table_data = []
        for sample, trimming_output in trimming_output_by_sample.items():
            # Add stats
            stats = trimming_output.informs_trimmomatic
            row = [sample.name_full, stats['paired_reads_in'], stats['paired_reads_out'], stats['forward_only_reads'],
                   stats['reverse_only_reads'], stats['reads_drop']]

            # Add FastQC reports
            for i, f in enumerate(trimming_output.fastq_reports_pre, start=1):
                relative_path = os.path.join('fastqc_report', f'{sample.name_valid}_{i}.html')
                section.add_file(f.path, relative_path)
                row.append(HtmlTableCell('view', link=relative_path))
            table_data.append(row)

        # Add table
        header = ['Sample', 'Total read pairs', 'Both surviving', 'Forward only surviving', 'Reverse only surviving',
                  'Dropped', 'FastQC report (fwd.)', 'FastQC report (rev.)']
        section.add_table(table_data, header, [('class', 'data')])

        # Save the report
        report.add_html_object(section)
        section.copy_files(report.output_dir)
        report.save()

    @staticmethod
    def construct_snp_matrix(sample_names: List[str], vcf_files: List[ToolIOFile], working_dir: str,
                             include_ref: bool = False) -> ToolIOFile:
        """
        Constructs a SNP matrix based on the given VCF files.
        :param sample_names: Sample names
        :param vcf_files: VCF files
        :param working_dir: Working directory
        :param include_ref: If True, the reference is included in the SNP matrix
        :return: SNP matrix ToolIOFile
        """
        if not os.path.isdir(working_dir):
            os.makedirs(working_dir)
        snp_matrix_constructor = SnpMatrixConstructor(Camel.get_instance())
        snp_matrix_constructor.update_parameters(include_ref='true' if include_ref else 'false')
        snp_matrix_constructor.add_input_files({
            'VCF': vcf_files,
            'SAMPLE_NAME': [ToolIOValue(s) for s in sample_names],
        })
        snp_matrix_constructor.run(working_dir)
        return snp_matrix_constructor.tool_outputs['FASTA'][0]

    @staticmethod
    def add_output_files_section(report: HtmlReport, column_names: List[str],
                                 output_files: Dict['SnpPhylogenyUtils.Sample', List[Union[str, None]]],
                                 snp_matrix: str):
        """
        Adds the section with the output files.
        :param report: Report
        :param column_names: Output file column names
        :param output_files: Output files
        :param snp_matrix: Path to the SNP matrix
        :return: None
        """
        section = HtmlReportSection('Output files')

        # Add SNP matrix
        relative_path = 'snp_matrix.fasta'
        section.add_file(snp_matrix, relative_path)
        table_data_snp_matrix = [
            ['SNP matrix:', HtmlTableCell('Download (FASTA)', link=relative_path)],
            ['Size:', SnpPhylogenyUtils.get_snp_matrix_size(snp_matrix)]
        ]
        section.add_table(table_data_snp_matrix, table_attributes=[('class', 'information')])

        # Add other output files table
        table_data = []
        for sample, output in output_files.items():
            row = [sample.name_full]
            for file_ in output:
                if file_ is None:
                    row.append('Not available')
                else:
                    relative_path = os.path.join(sample.name_valid, os.path.basename(file_))
                    row.append(HtmlTableCell('Download', link=section.add_file(file_, relative_path)))
                table_data.append(row)
        header = ['Sample'] + column_names
        section.add_table(table_data, header, [('class', 'data')])
        report.add_html_object(section)
        section.copy_files(report.output_dir)
        report.save()

    @staticmethod
    def add_metrics_section(report: HtmlReport, stats: Dict[Sample, List], header: List[str]) -> None:
        """
        Adds a section with analysis metrics to the report.
        :param report: Report
        :param stats: Stats
        :param header: Table header
        :return: None
        """
        section = HtmlReportSection('Analysis metrics')
        table_data = [[sample.name_full] + sts for sample, sts in stats.items()]
        section.add_table(table_data, header, [('class', 'data')])
        report.add_html_object(section)
        report.save()

    @staticmethod
    def get_snp_matrix_size(snp_matrix: str) -> int:
        """
        Returns the length of the given SNP matrix.
        Note: this function only checks the length of the first entry.
        :param snp_matrix: SNP matrix
        :return: SNP matrix size
        """
        with open(snp_matrix) as handle:
            entries = list(SeqIO.parse(handle, 'fasta'))
        return len(entries[0].seq)

    @staticmethod
    def check_snp_matrix_size(snp_matrix: str, size_min: int = 10, size_max: int = 25000) -> None:
        """
        Checks if the size of the SNP matrix is OK to perform model selection / tree building.
        :param snp_matrix: SNP matrix
        :param size_min: Minimum SNP matrix size
        :param size_max: Maximal SNP matrix size
        :return: None
        """
        snp_matrix_size = SnpPhylogenyUtils.get_snp_matrix_size(snp_matrix)
        if snp_matrix_size > size_max:
            raise ValueError('SNP matrix is too big ({}, max={}) to perform model selection / tree building'.format(
                snp_matrix_size, size_max))
        elif snp_matrix_size < size_min:
            raise ValueError('SNP matrix is too small ({}, min={}) to perform model selection / tree building'.format(
                    snp_matrix_size, size_min))

    @staticmethod
    def add_model_selection_section(report: HtmlReport, model_selection: Optional[ModelSelection] = None,
                                    error_message: Optional[str] = None) -> None:
        """
        Adds the section with the model selection results.
        :param report: Report
        :param model_selection: Model selection tool
        :param error_message: Error message
        :return: None
        """
        section = HtmlReportSection('Model selection')
        if error_message is not None:
            section.add_error_message(error_message)
        else:
            relative_path = 'model_selection_out.csv'
            section.add_file(model_selection.tool_outputs['CSV'][0].path, relative_path)
            table_data = [
                ['Selected model:', model_selection.informs['model_full']],
                ['Rates among sites:', model_selection.informs['rates_among_sites_full']],
                ['Overview:', HtmlTableCell('Download (CSV)', link=relative_path)]
            ]
            section.add_table(table_data, table_attributes=[('class', 'information')])
        report.add_html_object(section)
        section.copy_files(report.output_dir)
        report.save()

    @staticmethod
    def run_model_selection(snp_matrix: ToolIOFile, args: argparse.Namespace) -> ModelSelection:
        """
        Runs the MEGA model selection step.
        :param snp_matrix: SNP matrix
        :param args: Command line arguments
        :return: Model selection
        """
        model_selection = ModelSelection(Camel.get_instance())
        model_selection.add_input_files({'FASTA': [snp_matrix]})
        MEGAUtils.update_model_selection_parameters(
            model_selection, args.missing_data, args.branch_swap, args.site_cov_cutoff, args.threads)
        working_dir = os.path.join(args.working_dir, 'model_selection')
        if not os.path.isdir(working_dir):
            os.mkdir(working_dir)
        model_selection.run(working_dir)
        return model_selection

    @staticmethod
    def add_tree_building_section(report: HtmlReport, newick_path: Optional[str] = None,
                                  error_message: Optional[str] = None) -> None:
        """
        Adds the tree building section.
        :param report: Report
        :param newick_path: Path to the newick tree
        :param error_message: Error message
        :return: None
        """
        section = HtmlReportSection('Tree construction')
        if error_message is not None:
            section.add_error_message(error_message)
        else:
            # Add download link
            relative_path = 'tree.nwk'
            section.add_file(newick_path, relative_path)
            table_data = [['Tree:', HtmlTableCell('Download (Newick)', link=relative_path)]]
            section.add_table(table_data, table_attributes=[('class', 'information')])

            # Render tree
            output_path = os.path.join(report.output_dir, 'tree.png')
            NewickUtils.render(Camel(logging_config=None), newick_path, output_path, 'clad')
            section.add_html_object(HtmlElement('img', attributes=[('src', 'tree.png'), ('border', '1')]))

        report.add_html_object(section)
        section.copy_files(report.output_dir)
        report.save()

    @staticmethod
    def run_tree_building(snp_matrix: ToolIOFile, model: str, rates: str, args: argparse.Namespace) -> \
            MLTreeConstruction:
        """
        Builds a phylogenetic tree based on the given SNP matrix.
        :param snp_matrix: SNP matrix
        :param model: Selected model (e.g. 'T92')
        :param rates: Rates among sites (e.g. 'U')
        :param args: Command line arguments
        :return: Tree building tool instance
        """
        tree_building = MLTreeConstruction(Camel(logging_config=None))
        tree_building.add_input_files({'FASTA': [snp_matrix]})
        MEGAUtils.update_tree_building_parameters(
            tree_building, model, rates, args.bootstraps, args.missing_data, args.site_cov_cutoff,
            args.ml_method, args.branch_swap, args.threads)
        working_dir = os.path.join(args.working_dir, 'tree_building')
        if not os.path.isdir(working_dir):
            os.mkdir(working_dir)
        tree_building.run(working_dir)
        return tree_building
