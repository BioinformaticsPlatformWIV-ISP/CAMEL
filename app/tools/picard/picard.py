import logging
import os
import re

from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class Picard(Tool):
    """
    Super class for Picard tools
    """

    def __init__(self, tool_name, version, camel):
        """
        Initialize a picard tool
        :param tool_name: tool name
        :param version: Tool version
        :param camel: Camel instance
        :return: None
        """
        super(Picard, self).__init__(tool_name, version, camel)

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

    def _execute_tool(self):
        """
        Function to run Picard function
        :return: None
        """
        self._set_output()
        self._build_command()
        self._execute_command()
        self._set_inform()

    def _check_input(self):
        """
        Set the input specification, this default function handles only one SAM or BAM file as input
        :return: None
        """
        super(Picard, self)._check_input()

        self._check_required_inputs()

        if len(self._supported_inputs) != 0:
            input_type, input_files = self._check_supported_input()
            if len(input_files) > 1:
                raise ValueError("Can only specify one file of type {!r} as input of Picard {!r}.".format(
                    input_type, self._function_name))
            self._input_string = " I={}".format(input_files[0].path)

        self._set_input()

    def _check_required_inputs(self):
        """
        Check input requirements to run Picard function
        :return: None
        """
        for input_file in self._required_inputs:
            if input_file not in self._tool_inputs:
                raise KeyError('Picard {!r} required input file of type {!r} is missing!'.format(self._function_name, input_file))

    def _set_input(self):
        """
        Function to set required and optional inputs in self._input_string
        :return: None
        """
        if 'FASTA_REF' in self._tool_inputs:
            self._input_string += " R={}".format(self._tool_inputs['FASTA_REF'][0].path)

    def _check_supported_input(self):
        """
        Check supported (alternatives but still required) input to run Picard function
        :return: None
        """
        for input_type in self._supported_inputs:
            if input_type in self._tool_inputs:
                return input_type, self._tool_inputs[input_type]

        raise KeyError('None of the supported input types {!r} of Picard {!r} is specified!'.format(
            self._supported_inputs, self._function_name))

    def _set_output(self):
        """
        Set the output specification, this default function handles only one BAM file as output
        :return: None
        """
        self._tool_outputs['BAM'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output'].value))]

    def _build_command(self):
        """
        Build the command to run tool
        :return: None
        """
        self._command.command = " ".join([
            self._tool_command, self._input_string, self._output_string,
            " ".join(self._build_options(excluded_parameters=self._specific_parameters, delimiter='='))
        ])

    def _set_inform(self):
        """
        Analyse the result of picard run and update tool.informs, implement when necessary
        :return: None
        """
        pass

    def _check_command_output(self):
        """
        Analyse the result of Picard run
        :return: stdout_lines, standard output separated into lines
        """
        stdout_lines = self.stdout.split('\n')

        run_status = stdout_lines[-2].rstrip()
        if not re.match('Exit status: 0', run_status):
            raise ToolExecutionError("Picard {!r} fails to run, error msg: \n{}".format(
                self._function_name, self.stdout))

        # log WARNINGs in info.log
        for l in stdout_lines:
            if re.match('WARNING', l):
                logging.warning(" Picard - {}".format(l))

        return stdout_lines
