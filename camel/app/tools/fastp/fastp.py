import json
from pathlib import Path
from typing import Optional

from camel.app.camel import Camel
from camel.app.tools.tool import Tool
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile


class Fastp(Tool):
    """
    Fastp is an ultra-fast all-in-one FASTQ preprocessor.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('fastp', '0.23.4', camel)

    def _check_input(self) -> None:
        """
        Checks whether the required input files are specified.
        :return: None
        """
        if 'FASTQ' not in self._tool_inputs:
            raise InvalidInputSpecificationError(f'FASTQ input is required')
        super()._check_input()

    def __get_output_path(self, orientation: Optional[str], compress: bool = True) -> Path:
        """
        Returns the output filename.
        :param orientation: Read orientation
        :param compress: If True, output file is compressed
        :return: Path to output file
        """
        if orientation is None:
            return self.folder / (f"{self._parameters['output_name'].value}.fastq" + ('.gz' if compress else ''))
        return self.folder / (f"{self._parameters['output_name'].value}_{orientation}.fastq" + (
            '.gz' if compress else ''))

    def _build_command(self) -> None:
        """
        Build the command to run tool
        :return: None
        """
        parts = [self._tool_command]
        if len(self._tool_inputs['FASTQ']) == 1:
            parts.append(f"--in1 {self._tool_inputs['FASTQ'][0].path}")
            parts.append(f"--out1 {self.__get_output_path(orientation=None)}")
        elif len(self._tool_inputs['FASTQ']) == 2:
            parts.append(f"--in1 {self._tool_inputs['FASTQ'][0].path}")
            parts.append(f"--out1 {self.__get_output_path(orientation='1P')}")
            parts.append(f"--in2 {self._tool_inputs['FASTQ'][1].path}")
            parts.append(f"--out2 {self.__get_output_path(orientation='2P')}")
            parts.append(f"--unpaired1 {self.__get_output_path(orientation='1U')}")
            parts.append(f"--unpaired2 {self.__get_output_path(orientation='2U')}")
        else:
            raise InvalidInputSpecificationError(f'Invalid number of FASTQ input files (1 or 2 expected)')
        self._command.command = ' '.join([*parts, *self._build_options(excluded_parameters=['output_name'])])

    def _check_command_output(self) -> None:
        """
        Checks command output
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Tool execution failed: {self._command.stderr}")

    def _set_output(self) -> None:
        """
        set the output file to check
        """
        if len(self._tool_inputs['FASTQ']) == 1:
            self._tool_outputs['FASTQ'] = [ToolIOFile(self.__get_output_path(orientation=None))]
        else:
            self._tool_outputs['FASTQ_PE'] = [
                ToolIOFile(self.__get_output_path(orientation='1P')),
                ToolIOFile(self.__get_output_path(orientation='2P')),
            ]
            # Add SE reads
            if self.__get_output_path(orientation='1U').exists():
                self._tool_outputs['FASTQ_SE_FWD'] = [ToolIOFile(self.__get_output_path(orientation='1U'))]
            if self.__get_output_path(orientation='2U').exists():
                self._tool_outputs['FASTQ_SE_REV'] = [ToolIOFile(self.__get_output_path(orientation='2U'))]

        # General output files (these are always created)
        self._tool_outputs['JSON'] = [ToolIOFile(self.folder / 'fastp.json')]
        self._tool_outputs['HTML'] = [ToolIOFile(self.folder / 'fastp.html')]

    def _set_informs(self) -> None:
        """
        Sets the tool informs.
        :return: None
        """
        with self._tool_outputs['JSON'][0].path.open() as handle:
            data = json.load(handle)
        self._informs['summary'] = data['summary']

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._set_output()
        self._set_informs()
