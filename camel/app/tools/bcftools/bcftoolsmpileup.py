from pathlib import Path
from typing import Optional

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftoolsbase import BcftoolsBase
from camel.app.tools.toolpipeable import ToolPipeable


class BcftoolsMpileup(BcftoolsBase, ToolPipeable):
    """
    Multi-way pileup producing genotype likelihoods.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('bcftools mpileup', '1.17', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Reference genome input is required (FASTA)")
        if 'BAM' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Alignment input is required (BAM)")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command(self._get_output_path())
        self._execute_command()
        self._tool_outputs[self._get_output_key()] = [ToolIOFile(self._get_output_path())]

    def _build_command(self, path_out: Optional[Path], pipe_in: bool = False, pipe_out: bool = False) -> None:
        """
        Builds the command that is called.
        :param pipe_in: True if the tool is the input of a pipe
        :param pipe_out: True if the tool is the output of a pipe
        :param path_out: Output path
        :return: None
        """
        command_parts = [
            self._tool_command,
            str(self._tool_inputs['BAM'][0].path) if not pipe_in else '-',
            f"--fasta-ref {self._tool_inputs['FASTA'][0].path}"
        ]
        if 'BED_include' in self._tool_inputs:
            command_parts.append(f"--targets-file {self._tool_inputs['BED_include'][0].path}")
        elif 'BED_exclude' in self._tool_inputs:
            command_parts.append(f"--targets-file ^{self._tool_inputs['BED_exclude'][0].path}")
        command_parts.extend(self._build_options(excluded_parameters=['output_filename']))
        if not pipe_out:
            command_parts.extend(['--output', str(path_out)])
        self._command.command = ' '.join(command_parts)

    def _before_pipe(self, dir_, pipe_in: bool, pipe_out: bool) -> None:
        """
        Prepares the command that will be piped.
        :param dir_: Running directory
        :param pipe_in: True if tool receives piped input
        :param pipe_out: True if tool generates piped output
        :return: None
        """
        self._build_command(self._get_output_path(), pipe_in, pipe_out)

    def _after_pipe(self, stderr: str, is_last_in_pipe: bool) -> None:
        """
        Performs the required steps after executing the tool as part of a pipe.
        :param stderr: Stderr for this command in the pipe
        :param is_last_in_pipe: Boolean to indicate if this is the last step in the pipe
        :return: None
        """
        if is_last_in_pipe:
            self._tool_outputs[self._get_output_key()] = [ToolIOFile(self._get_output_path())]
