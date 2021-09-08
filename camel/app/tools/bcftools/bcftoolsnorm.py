import re
from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftools import Bcftools
from camel.app.tools.samtools.samtools import Samtools


class BcftoolsNorm(Bcftools):
    """
    Indexes bgzip compressed VCF files and BCF files.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('bcftools norm', '1.9', camel)
        self._input_key = None

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('VCF', 'VCF_GZ')):
            raise InvalidInputSpecificationError("No input file found (VCF / VCF_GZ supported)")
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Reference genome input is required (FASTA).")
        super(BcftoolsNorm, self)._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        output_path = self._get_output_path()
        self.__build_command(output_path)
        self._execute_command()
        self.__set_output(output_path)
        self.__set_informs()

    def __build_command(self, output_path: Path) -> None:
        """
        Builds the command for this tool.
        :param output_path: Path to the output file
        :return: None
        """
        input_key = next(key for key in ('VCF', 'VCF_GZ') if key in self._tool_inputs)
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(self._build_options()),
            f"-f {self._tool_inputs['FASTA'][0].path}",
            str(self._tool_inputs[input_key][0].path),
            f' > {output_path}'
        ])

    def __set_output(self, output_path: Path) -> None:
        """
        Sets the output of this tool.
        :param output_path: Output path
        :return: None
        """
        output_key = self._get_output_key()
        self._tool_outputs[output_key] = [ToolIOFile(output_path)]

    def _check_command_output(self) -> None:
        """
        Checks if the tool was executed successfully.
        :return: None
        """
        if not self._command.returncode == 0:
            raise ToolExecutionError(self._command.stderr)

    def __set_informs(self) -> None:
        """
        Sets the informs for this tool.
        :return: None
        """
        for line in self.stderr.splitlines():
            m = re.search(r'total/split/realigned/skipped:\t(\d+)/(\d+)/(\d+)/(\d+)', line)
            if m is None:
                continue
            self._informs['total'] = int(m.group(1))
            self._informs['split'] = int(m.group(2))
            self._informs['realigned'] = int(m.group(3))
            self._informs['skipped'] = int(m.group(4))
            return
        raise ToolExecutionError("Cannot extract informs")
