import shutil

import logging
import traceback

import os

import abc
import argparse

from Bio import SeqIO

from app.camel import Camel
from app.components.files.fastqutils import FastqUtils
from app.components.vcf import vcfutils
from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliodirectory import ToolIODirectory
from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue
from app.pipeline.pipeline import Pipeline
from app.tools.mega.mltreeconstruction import MLTreeConstruction
from app.tools.mega.modelselection import ModelSelection
from app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
from app.tools.treevector.treevector import TreeVector
from resources import YAML_READ_TRIMMING
from scripts.snpphylogeny.htmlreportersnpphylogeny import HtmlReporterSnpPhylogeny


class SnpPhylogeny(object):
    """
    Base class for the SNP phylogeny pipeline scripts.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, name):
        """
        Initializes the pipeline script.
        :param name: Pipeline name
        """
        self._name = name
        self._camel = Camel()
        self._destination_path = os.path.abspath('.')
        self._args = self._parse_arguments()
        self._html = HtmlReporterSnpPhylogeny(self._name, self._args.dir_html)
        self._sample_names = self.__get_sample_names()
        self._bam_files = []
        self._vcf_files = []
        self._snp_pipeline_input = self.__get_input_reads()

    def _parse_arguments(self):
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--html', required=True)
        argument_parser.add_argument('--dir-html', required=True)
        argument_parser.add_argument('--reference', required=True)
        argument_parser.add_argument('--reference-name')
        argument_parser.add_argument('--sample', nargs=5, action='append', required=True)
        argument_parser.add_argument('--trim-reads', action='store_true')
        argument_parser.add_argument('--missing-data', choices=['complete_deletion', 'all_sites', 'partial_deletion'])
        argument_parser.add_argument('--site-cov-cutoff', choices=range(0, 101), type=int)
        argument_parser.add_argument('--branch-swap', choices=['none', 'weak', 'very_weak', 'moderate', 'strong',
                                                               'very_strong'])
        argument_parser.add_argument('--bootstraps', type=int, required=True)
        argument_parser.add_argument('--ml-method', choices=['nni', 'spr3', 'spr5'], required=True)
        self._add_specific_arguments(argument_parser)
        return argument_parser.parse_args()

    def __check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if len(self._args.sample) < 4:
            raise RuntimeError("At least 4 samples are required ")

    @staticmethod
    def _add_specific_arguments(argument_parser):
        """
        Parses the pipeline specific arguments, should be implemented by the subclasses.
        :return: None
        """
        pass

    def run(self):
        """
        Runs the SNP phylogeny pipeline script.
        :return: None
        """
        try:
            os.mkdir(self._args.dir_html)
            self.__check_input()
            # Initialize HTML report
            self._html.initialize()
            with open(self._args.html, 'w') as handle:
                handle.write(self._html.get_html())
            # Top sections
            self._html.add_analysis_info_section(len(self._args.sample), self._args.reference_name)
            self._html.add_parameter_section(self._args)
            # Read trimming
            if self._args.trim_reads:
                self.__run_read_trimming()
            # SNP Calling
            snp_matrix, output_files = self._run_snp_calling(self._snp_pipeline_input)
            # Metrics
            metrics = self._collect_metrics()
            self._html.add_metrics_section(metrics)
            try:
                size = self.__get_snp_matrix_size(snp_matrix)
                self._html.add_output_files_section(self._sample_names, snp_matrix, size, output_files)
                if size < 5:
                    raise ValueError(
                        "SNP matrix is too small ({} positions) to construct a phylogenetic tree.".format(size))
                # Tree building
                model, rates = self._run_model_selection(snp_matrix)
                try:
                    self._run_tree_building(snp_matrix, model, rates)
                except ToolExecutionError:
                    self._html.add_error_message(
                        """Could not build bootstrap tree, check the logs for more details. 
                        The SNP matrix might be too small, try using less stringent filters.""")
            except ValueError as err:
                logging.info(traceback.format_exc())
                self._html.add_error_message(err.message)
            # Save report
            with open(self._args.html, 'w') as handle:
                handle.write(self._html.get_html())
        except:
            raise
        finally:
            self.__copy_log_file()

    def __get_sample_names(self):
        """
        Retrieves the sample names from the command line arguments. This function also updates the sample names in 
        self._args. 
        :return: None
        """
        sample_names = []
        for i in range(0, len(self._args.sample)):
            if self._args.sample[i][0] == '':
                sample_name = FastqUtils.get_sample_name(os.path.basename(self._args.sample[i][1]))
                sample_names.append(sample_name)
                self._args.sample[i][0] = sample_name
            else:
                sample_names.append(self._args.sample[i][0])
        return sample_names

    def __get_input_reads(self):
        """
        Returns the default input for the SNP pipeline.
        :return: Input reads dictionary
        """
        return [{'FASTQ_PE': [ToolIOFile(fwd_data), ToolIOFile(rev_data)]} for
                _, _, fwd_data, _, rev_data in self._args.sample]

    def __run_read_trimming(self):
        """
        Runs the read trimming pipeline on all samples.
        :return: None
        """
        trimming_pipelines = {}
        for (sample_name, _, forward_data, _, reverse_data), i in zip(
                self._args.sample, range(0, len(self._args.sample))):
            pipeline = Pipeline([YAML_READ_TRIMMING], self._camel, True)
            pipeline.set_initial_input({
                'FASTQ': [ToolIOFile(forward_data), ToolIOFile(reverse_data)],
                'SAMPLE_NAME': [ToolIOValue(sample_name)],
                'HTML': [ToolIOFile(self._args.html)],
                'DIR_HTML': [ToolIODirectory(self._args.dir_html)]
            })
            pipeline.set_configs({'SKIP_REPORT': True})
            pipeline.add_job_options({'Report_generation': {'disable_html_output': True}})
            pipeline.run(self._destination_path)
            trimming_pipelines[sample_name] = pipeline
            self._snp_pipeline_input[i] = pipeline.get_step('Read_trimming').outputs
        self._html.add_trimming_section(trimming_pipelines)

    @abc.abstractmethod
    def _run_snp_calling(self, reads):
        """
        Runs the SNP calling pipeline. Has to be implemented by the subclasses.
        :param reads: PE reads for all samples
        :return: SNP matrix, output files
        :rtype: ToolIOFile, [ToolIOFile]
        """
        return

    @staticmethod
    def __get_snp_matrix_size(snp_matrix):
        """
        Returns the length of the given SNP matrix.
        :param snp_matrix: SNP matrix
        :return: SNP matrix size
        """
        with open(snp_matrix.path, 'r') as handle:
            entries = list(SeqIO.parse(handle, 'fasta'))
        return len(entries[0].seq)

    def _collect_metrics(self):
        """
        Collects analysis statistics.
        :return: None
        """
        metrics = [['Sample', 'Total reads', 'Mapping rate', 'SNPs']]
        for sample_name, bam_file, vcf_file in zip(self._sample_names, self._bam_files, self._vcf_files):
            samtools_flagstat = SamtoolsFlagstat(self._camel)
            samtools_flagstat.add_input_files({'BAM': [bam_file]})
            samtools_flagstat.run(self._destination_path)
            total_reads = int(samtools_flagstat.informs['total'][0])
            mapping_rate = float(samtools_flagstat.informs['mapped'][0]) / total_reads
            snps = vcfutils.count_variants(vcf_file.path)
            metrics.append([sample_name, total_reads, '{:.2f}'.format(mapping_rate*100), snps])
        return metrics

    def _run_model_selection(self, snp_matrix):
        """
        Runs the MEGA model selection step.
        :param snp_matrix: SNP matrix
        :return: None
        """
        model_selection = ModelSelection(self._camel)
        model_selection.add_input_files({'FASTA': [snp_matrix]})
        if self._args.missing_data == 'complete_deletion':
            model_selection.update_parameters(missing_data_treatment='Complete deletion')
        elif self._args.missing_data == 'all_sites':
            model_selection.update_parameters(missing_data_treatment='Use all sites')
        elif self._args.missing_data == 'partial_deletion':
            model_selection.update_parameters(missing_data_treatment='Partial deletion',
                                              site_coverage_cutoff=self._args.site_cov_cutoff)
        model_selection.update_parameters(branch_swap_filter=self._args.branch_swap.title().replace('_', ' '))
        model_selection_dir = os.path.join(self._destination_path, 'model_selection')
        os.mkdir(model_selection_dir)
        model_selection.run(model_selection_dir)
        self._html.add_model_selection_section(model_selection)
        return model_selection.informs['model'], model_selection.informs['rates_among_sites']

    def _run_tree_building(self, snp_matrix, model, rates):
        """
        Builds a phylogenetic tree based on the given SNP matrix.
        :param snp_matrix: SNP matrix
        :param model: Selected model
        :param rates: Rates among sites
        :return: Tree building instance
        """
        tree_building = MLTreeConstruction(self._camel)
        tree_building.add_input_files({'FASTA': [snp_matrix]})
        tree_building.update_parameters(bootstrap_replications=self._args.bootstraps,
                                        test_of_phylogeny='Bootstrap method')
        if rates == 'G+I':
            tree_building.update_parameters(rates_among_sites='G+I')
            tree_building.update_parameters(gamma_categories='5')
        elif rates == 'G':
            tree_building.update_parameters(rates_among_sites='G')
            tree_building.update_parameters(gamma_categories='5')
        elif rates == 'I':
            tree_building.update_parameters(rates_among_sites='I')
        else:
            tree_building.update_parameters(rates_among_sites='U')

        if self._args.missing_data == 'complete_deletion':
            tree_building.update_parameters(missing_data_treatment='Complete deletion')
        elif self._args.missing_data == 'all_sites':
            tree_building.update_parameters(missing_data_treatment='Use all sites')
        elif self._args.missing_data == 'partial_deletion':
            tree_building.update_parameters(missing_data_treatment='Partial deletion',
                                            site_coverage_cutoff=self._args.site_cov_cutoff)

        tree_building.update_parameters(model=model)
        tree_building.update_parameters(branch_swap_filter=self._args.branch_swap.title().replace('_', ' '))
        tree_building.update_parameters(heuristic_method=self._args.ml_method.upper())
        tree_building.run(self._destination_path)
        newick_tree = self.__render_tree(tree_building.tool_outputs['NWK'][0])
        self._html.add_tree_building_section(tree_building, newick_tree, self._args.bootstraps)
        return tree_building

    def __render_tree(self, newick_tree):
        """
        Renders the Newick tree to image format.
        :param newick_tree: Newick tree
        :return: Rendered tree
        """
        tree_vector = TreeVector(self._camel)
        tree_vector.add_input_files({'NWK': [newick_tree]})
        tree_vector.update_parameters(output_format='png', output_filename='tree.png', type='clad')
        tree_vector.run()
        return tree_vector.tool_outputs['PNG'][0]

    def __copy_log_file(self):
        """
        Tries to copy the CAMEL log file to the output directory.
        :return: None
        """
        log_path = os.path.join(self._destination_path, 'camel.log')
        if os.path.isfile(log_path):
            shutil.copy(log_path, os.path.join(self._args.dir_html, 'log.txt'))
