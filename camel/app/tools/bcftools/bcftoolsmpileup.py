from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftoolsbase import BcftoolsBase


class BcftoolsMpileup(BcftoolsBase):
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

    def _build_command(self, path_out: Path) -> None:
        """
        Builds the command that is called.
        :param path_out: Output path
        :return: None
        """
        command_parts = [
            self._tool_command,
            str(self._tool_inputs['BAM'][0].path),
            f"--fasta-ref {self._tool_inputs['FASTA'][0].path}"
        ]
        if 'BED_include' in self._tool_inputs:
            command_parts.append(f"--targets-file {self._tool_inputs['BED_include'][0].path}")
        elif 'BED_exclude' in self._tool_inputs:
            command_parts.append(f"--targets-file ^{self._tool_inputs['BED_exclude'][0].path}")
        command_parts.extend(self._build_options(excluded_parameters=['output_filename']))
        command_parts.extend(['--output', str(path_out)])
        self._command.command = ' '.join(command_parts)
