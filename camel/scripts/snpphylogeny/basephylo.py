import argparse
import logging
import os
from typing import Optional, List, Dict

import abc

from camel.app.components.phylogeny.snpphylogenyutils import SnpPhylogenyUtils, InvalidInputError
from camel.app.components.workflows.readtrimmingwrapper import ReadTrimmingWrapper
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.mega.modelselection import ModelSelection


class BasePhylo(object, metaclass=abc.ABCMeta):
    """
    Base class for the SNP phylogeny pipelines.
    """

    def __init__(self, pipeline_name: str, args: Optional[argparse.Namespace]) -> None:
        """
        Initializes the base class.
        :param pipeline_name: Pipeline name (e.g. 'Samtools')
        :param args: Command line arguments
        """
        self._pipeline_name = pipeline_name
        self._args = args if args is not None else self._parse_arguments()
        self._report = SnpPhylogenyUtils.initialize_report(self._pipeline_name, self._args)
        self._samples = self.__extract_samples()

    @property
    def samples_by_name(self) -> Dict[str, SnpPhylogenyUtils.Sample]:
        """
        Returns the samples as a dictionary with the sample name as key.
        :return: Samples by name
        """
        return {s.name_valid: s for s in self._samples}

    @staticmethod
    @abc.abstractmethod
    def _parse_arguments() -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        pass

    def __extract_samples(self) -> List[SnpPhylogenyUtils.Sample]:
        """
        Extracts sample objects from the provided input.
        :return: List of samples
        """
        try:
            return SnpPhylogenyUtils.extract_samples(self._args)
        except InvalidInputError as err:
            logging.error(f"Invalid input: {err}")
            self._report.add_error_message(str(err))
            self._report.save()
            exit(0)

    def _get_mapping_input(self) -> Dict[SnpPhylogenyUtils.Sample, SnpPhylogenyUtils.MappingInput]:
        """
        Returns the input for the read mapping.
        :return: Mapping input per sample
        """
        if self._args.trim_reads:
            mapping_input_by_sample = {}
            trimming_output_by_sample = {}
            for sample in self._samples:
                working_dir = os.path.join(self._args.working_dir, sample.name_valid, 'trimming')
                wrapper = ReadTrimmingWrapper(working_dir)
                wrapper.run_workflow([f.path for f in sample.reads_raw])
                # noinspection PyCallByClass
                mapping_input_by_sample[sample] = SnpPhylogenyUtils.MappingInput(
                    pe=wrapper.output.trimmed_reads_pe,
                    se_fwd=wrapper.output.trimmed_reads_se_fwd[0],
                    se_rev=wrapper.output.trimmed_reads_se_rev[0]
                )
                trimming_output_by_sample[sample] = wrapper.output
            SnpPhylogenyUtils.add_trimming_section(self._report, trimming_output_by_sample)
            return mapping_input_by_sample
        else:
            SnpPhylogenyUtils.add_trimming_section_empty(self._report)
            # noinspection PyCallByClass
            return {s: SnpPhylogenyUtils.MappingInput(pe=s.reads_raw) for s in self._samples}

    def _run_model_selection(self, snp_matrix: ToolIOFile) -> ModelSelection:
        """
        Runs the model selection.
        Quits the programs if the SNP matrix is too small.
        :param snp_matrix: SNP matrix
        :return: Model selection tool instance
        """
        try:
            SnpPhylogenyUtils.check_snp_matrix_size(snp_matrix.path)
        except ValueError as err:
            SnpPhylogenyUtils.add_model_selection_section(self._report, error_message=str(err))
            exit(0)
        else:
            model_selection = SnpPhylogenyUtils.run_model_selection(snp_matrix, self._args)
            SnpPhylogenyUtils.add_model_selection_section(self._report, model_selection=model_selection)
            return model_selection

    def _run_tree_building(self, snp_matrix: ToolIOFile, model_selection: ModelSelection) -> None:
        """
        Runs the tree building.
        :param snp_matrix: SNP matrix
        :param model_selection: Model selection instance.
        :return: None
        """
        try:
            tree_building = SnpPhylogenyUtils.run_tree_building(
                snp_matrix, model_selection.informs['model'], model_selection.informs['rates_among_sites'], self._args)
            SnpPhylogenyUtils.add_tree_building_section(self._report, tree_building.tool_outputs['NWK'][0].path)
        except ToolExecutionError:
            SnpPhylogenyUtils.add_tree_building_section(
                self._report, error_message='Error constructing tree, SNP matrix might be too small')
