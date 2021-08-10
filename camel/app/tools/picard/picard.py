import abc
import logging
import re
from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Picard(Tool, metaclass=abc.ABCMeta):

    """
    Super class for Picard tools
    """

    def __init__(self, tool_name: str, version: str, camel: Camel) -> None:
        """
        Initialize a picard tool
        :param tool_name: tool name
        :param version: Tool version
        :param camel: Camel instance
        :return: None
        """
        super().__init__(tool_name, version, camel)

        self._function_name = None
        # parameters that should not be handled by self.build_options function
        self._specific_parameters = []
        self._required_inputs = ['BAM', 'SAM']
        self._main_input = ''
        self._input_string = ''
        self._output_string = ''
        self._output_type = 'BAM'
        # Elements for building command
        self._java_options = '-mx8G -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dpicard.useLegacyParser=false'
        self._java_options_temp_dir = 'TMP_DIR=/temp/picard'

    def update_java_options(self, java_options: str) -> None:
        """
        Returns the formatted java options of this tool.
        :return: Name
        """
        logging.info(f"Java options updated: '{java_options}'")
        self._java_options = f'{java_options}'

    def _execute_tool(self) -> None:
        """
        Function to run Picard function
        :return: None
        """
        self._set_input()
        self._set_output()
        self._build_command()
        self._execute_command()
        self._set_informs()

    def _check_input(self):
        """
        Check input for a tool and prepare command line parameters for input
        :return: None
        """
        # BAM or SAM input required; mutually exclusive
        if ('SAM' in self._tool_inputs) and ('BAM' in self._tool_inputs):
            raise InvalidInputSpecificationError()
        elif 'SAM' in self._tool_inputs:
            self._main_input = 'SAM'
            self._required_inputs.remove('BAM')
        elif 'BAM' in self._tool_inputs:
            self._main_input = 'BAM'
            self._required_inputs.remove('SAM')

        for input_file in self._required_inputs:
            if input_file not in self._tool_inputs:
                raise InvalidInputSpecificationError('Picard {!r} required {!r} input is missing in _tool_inputs!'.format(
                    self._name, input_file))

        super(Picard, self)._check_input()

    def _set_input(self) -> None:
        """
        Function to set main and extra inputs in self._input_string
        :return: None
        """
        if 'BAM' in self._tool_inputs:
            self._input_string += f"INPUT={self._tool_inputs['BAM'][0].path} "

        if 'SAM' in self._tool_inputs:
            self._input_string += f"INPUT={self._tool_inputs['SAM'][0].path} "

        if 'FASTA_REF' in self._tool_inputs:
            self._input_string += f"R={self._tool_inputs['FASTA_REF'][0].path} "

    def _set_output(self) -> None:
        """
        Set the output specification, this default function handles only one BAM file as output
        :return: None
        """
        self._tool_outputs[self._output_type] = [
            ToolIOFile(Path(self._folder) / self._parameters['output'].value)
        ]

    def _build_command(self) -> None:
        """
        Build the command to run tool
        :return: None
        """
        self._command.command = " ".join([
            "java", self._java_options, "-jar $PICARD_JAR", self._tool_command, self._java_options_temp_dir, self._input_string, self._output_string,
            ' '.join(self._build_options(excluded_parameters=self._specific_parameters, delimiter='=')), '2>&1'
        ])

    def _set_informs(self) -> None:
        """
        Analyse the result of picard run and update tool.informs, implement when necessary
        :return: None
        """
        pass

    def _check_command_output(self) -> None:
        """
        Verify tool execution (return code) and stdout
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f'Picard {self._function_name,} fails to run, error msg: \n{self.stdout}')

        # log WARNINGs
        for line in self.stdout.splitlines():
            if re.match('WARNING', line):
                logging.warning(f' Picard - {line}')
