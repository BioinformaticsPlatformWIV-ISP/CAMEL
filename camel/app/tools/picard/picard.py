import abc
import logging
import os
import re
from typing import Tuple, List

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
        # alternative types of files that can be used as main input of a picard tool
        # - reads mapping input: 'SAM', 'BAM'
        # - variance input: 'VCF' 'BCF'
        # ...
        self._supported_inputs = ['SAM', 'BAM']
        # individual files of different types that is required: e.g, FASTA_REF
        self._required_inputs = []
        self._input_string = ''
        self._output_string = ''

    def _execute_tool(self) -> None:
        """
        Function to run Picard function
        :return: None
        """
        self._set_output()
        self._build_command()
        self._execute_command()
        self._set_informs()

    def _check_input(self) -> None:
        """
        Set the input specification, this default function handles only one SAM or BAM file as input
        :return: None
        """
        super(Picard, self)._check_input()

        self._check_required_inputs()

        if len(self._supported_inputs) != 0:
            input_type, input_files = self._check_supported_input()
            if len(input_files) > 1:
                raise InvalidInputSpecificationError(f'Can only specify one file of type {input_type} '
                                                     f'as input of Picard {self._function_name}.')
            self._input_string = f' I={input_files[0].path}'

        self._set_input()

    def _check_required_inputs(self) -> None:
        """
        Check input requirements to run Picard function
        :return: None
        """
        for input_file in self._required_inputs:
            if input_file not in self._tool_inputs:
                raise InvalidInputSpecificationError(
                    f'Picard {self._function_name} required input file of type {input_file} is missing!')

    def _set_input(self) -> None:
        """
        Function to set required and optional inputs in self._input_string
        :return: None
        """
        if 'FASTA_REF' in self._tool_inputs:
            self._input_string += f" R={self._tool_inputs['FASTA_REF'][0].path}"

    def _check_supported_input(self) -> Tuple[str, List[ToolIOFile]]:
        """
        Check supported (alternatives but still required) input to run Picard function
        :return: None
        """
        for input_type in self._supported_inputs:
            if input_type in self._tool_inputs:
                return input_type, self._tool_inputs[input_type]

        raise InvalidInputSpecificationError(f'None of the supported input types {self._supported_inputs} '
                                             f'of Picard {self._function_name} is specified!')

    def _set_output(self) -> None:
        """
        Set the output specification, this default function handles only one BAM file as output
        :return: None
        """
        self._tool_outputs['BAM'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output'].value))]

    def _build_command(self) -> None:
        """
        Build the command to run tool
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command, self._input_string, self._output_string,
            ' '.join(self._build_options(excluded_parameters=self._specific_parameters, delimiter='='))
        ])

    def _set_informs(self) -> None:
        """
        Analyse the result of picard run and update tool.informs, implement when necessary
        :return: None
        """
        pass

    def _check_command_output(self) -> None:
        """
        Analyse the result of Picard run
        :return: None
        """
        stdout_lines = self.stdout.splitlines()

        run_status = stdout_lines[-1]
        if not re.match('Exit status: 0', run_status):
            raise ToolExecutionError(f'Picard {self._function_name,} fails to run, error msg: \n{self.stdout}')

        # log WARNINGs
        for line in stdout_lines:
            if re.match('WARNING', line):
                logging.warning(f' Picard - {line}')
