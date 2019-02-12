from typing import List

import os

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool
from camel.resources.r import RSCRIPT_QC_COVERAGE, RSCRIPT_QC_ST, RSCRIPT_QC_MAPPING_RATE


class QCVisualization(Tool):
    """
    This tool is used to visualize the results of the advanced QC checks.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        """
        super().__init__('QC Visualization', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'cov' not in self._input_informs:
            raise InvalidInputSpecificationError("Coverage informs required ('cov')")
        if 'mlst' not in self._input_informs:
            raise InvalidInputSpecificationError("(cg)MLST informs required ('mlst')")
        if 'map' not in self._input_informs:
            raise InvalidInputSpecificationError("Mapping informs required ('map')")
        super()._check_input()

    def __run_rscript(self, script: str, output_file: str, params: List[str]) -> None:
        """
        Runs R script to create the visualization.
        :param script: Script to run
        :param output_file: Output file
        :param params: Script parameters
        :return: None
        """
        self._command.command = f"Rscript --vanilla {script} {output_file} {' '.join(params)}"
        self._execute_command()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Coverage
        output_file = os.path.join(self._folder, 'coverage_plot.png')
        coverage = min(self._input_informs['cov']['median_depth'], 50)
        self.__run_rscript(RSCRIPT_QC_COVERAGE, output_file, ['{:.2f}'.format(coverage)])
        self._tool_outputs['PNG_cov'] = [ToolIOFile(output_file)]

        # cgMLST
        output_file = os.path.join(self._folder, 'typing.png')
        cgmlst_perc = '{:.2f}'.format(
            100 * self._input_informs['mlst']['hits_found'] / self._input_informs['mlst']['nb_of_loci'])
        self.__run_rscript(RSCRIPT_QC_ST, output_file, [cgmlst_perc, '"{}"'.format(
            self._input_informs['mlst']['title'])])
        self._tool_outputs['PNG_st'] = [ToolIOFile(output_file)]

        # Mapping rate
        output_file = os.path.join(self._folder, 'mapping_rate.png')
        mapping_rate = '{:.2f}'.format(float(self._input_informs['map']['stats_map_rate']))
        self.__run_rscript(RSCRIPT_QC_MAPPING_RATE, output_file, [mapping_rate])
        self._tool_outputs['PNG_mapping'] = [ToolIOFile(output_file)]

    def _check_command_output(self) -> None:
        """
        Checks if the tool execution was successful.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Error executing {self.name}:\n{self._command.stderr}")
