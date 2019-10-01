#!/usr/bin/env python
import argparse
from typing import Tuple, List, Optional

import os
import shutil

from camel.app.camel import Camel
from camel.app.components.phylogeny.megautils import MEGAUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.mega.mltreeconstruction import MLTreeConstruction
from camel.app.tools.mega.modelselection import ModelSelection
from camel.app.tools.snpmatrix.snpmatrixconstructor import SnpMatrixConstructor


class MainMega(object):
    """
    This class contains the main script for the MEGA model selection and tree building.
    """

    SNP_MATRIX_FILENAME = 'snp_matrix.fasta'

    def __init__(self, args: Optional[argparse.Namespace] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = args if args is not None else MainMega._parse_arguments()
        self._camel = Camel()

    @staticmethod
    def _parse_arguments() -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--fasta', help="Input SNP matrix FASTA file")
        argument_parser.add_argument('--vcf', action='append', nargs=2)
        argument_parser.add_argument('--output-tree', help="Output Newick tree file")
        argument_parser.add_argument('--output-model', help="Output tabular model selection file")
        argument_parser.add_argument('--output-snp-matrix',
                                     help='If set, the SNP matrix is exported to this FASTA file')
        argument_parser.add_argument('--action', choices=['both', 'model', 'tree'], required=True)
        argument_parser.add_argument('--missing-data', required=True,
                                     choices=['use_all_sites', 'complete_deletion', 'partial_deletion'])
        argument_parser.add_argument('--branch-swap', required=True, choices=['none', 'very_strong', 'strong',
                                                                              'moderate', 'weak', 'very_weak'])
        argument_parser.add_argument('--bootstraps', required=True, help="Number of bootstrap replications.",
                                     type=int)
        argument_parser.add_argument('--ml-method', required=True, choices=['nni', 'spr3', 'spr5'])
        argument_parser.add_argument('--site-cov-cutoff', choices=range(0, 101), type=int)
        argument_parser.add_argument('--model', choices=['JC', 'K2', 'T92', 'HKY', 'TN93', 'GTR'])
        argument_parser.add_argument('--rates', choices=['G+I', 'G', 'I', 'U'], default='U')
        argument_parser.add_argument('--working-dir', default=os.path.abspath('.'))
        argument_parser.add_argument('--threads', type=int, default=3)
        return argument_parser.parse_args()

    def run(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        # Prepare SNP matrix
        if self._args.fasta is not None:
            fasta_path = os.path.join(self._args.working_dir, MainMega.SNP_MATRIX_FILENAME)
            os.symlink(self._args.fasta, fasta_path)
        else:
            fasta_path = self.__build_snp_matrix(self._args.vcf)
            if self._args.output_snp_matrix is not None:
                shutil.copyfile(fasta_path, self._args.output_snp_matrix)

        # Perform required actions
        if self._args.action == 'model':
            self.__run_model_selection(fasta_path)
        elif self._args.action == 'tree':
            self.__run_tree_construction(fasta_path, self._args.model, self._args.rates)
        else:
            model, rates = self.__run_model_selection(fasta_path)
            self.__run_tree_construction(fasta_path, model, rates)

    def __build_snp_matrix(self, vcf_files: List[str]) -> str:
        """
        Builds a SNP matrix by combining all SNP positions across multiple VCF files.
        :param vcf_files: VCF files
        :return: SNP matrix path
        """
        snp_matrix_constructor = SnpMatrixConstructor(self._camel)
        snp_matrix_constructor.update_parameters(include_ref=None)
        snp_matrix_constructor.add_input_files({
            'VCF': [ToolIOFile(v) for v, _ in vcf_files],
            'SAMPLE_NAME': [ToolIOValue(n) for _, n in vcf_files]})
        snp_matrix_constructor.run(self._args.working_dir)
        return snp_matrix_constructor.tool_outputs['FASTA'][0].path

    def __run_model_selection(self, fasta_path: str) -> Tuple[str, str]:
        """
        Runs the model selection.
        :param fasta_path: FASTA input file
        :return: Model, rates parameter
        """
        model_selection = ModelSelection(self._camel)
        MEGAUtils.update_model_selection_parameters(
            model_selection, self._args.missing_data, self._args.branch_swap, self._args.site_cov_cutoff,
            self._args.threads)
        model_selection.add_input_files({'FASTA': [ToolIOFile(fasta_path)]})
        working_dir = os.path.join(self._args.working_dir, 'model_selection')
        if not os.path.isdir(working_dir):
            os.mkdir(working_dir)
        model_selection.run(working_dir)
        shutil.copyfile(model_selection.tool_outputs['CSV'][0].path, self._args.output_model)
        return model_selection.informs['model'], model_selection.informs['rates_among_sites']

    def __run_tree_construction(self, fasta_path: str, model: str, rates: str) -> None:
        """
        Runs the tree construction.
        :param fasta_path: FASTA input file
        :param model: Nucleotide substitution model
        :param rates: Rates among sites
        :return: None
        """
        tree_building = MLTreeConstruction(self._camel)
        tree_building.add_input_files({'FASTA': [ToolIOFile(fasta_path)]})
        MEGAUtils.update_tree_building_parameters(
            tree_building, model, rates, self._args.bootstraps, self._args.missing_data, self._args.site_cov_cutoff,
            self._args.ml_method, self._args.branch_swap, self._args.threads)
        working_dir = os.path.join(self._args.working_dir, 'tree_building')
        if not os.path.isdir(working_dir):
            os.mkdir(working_dir)
        tree_building.run(working_dir)
        shutil.copyfile(tree_building.tool_outputs['NWK'][0].path, self._args.output_tree)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainMega()
    main.run()
