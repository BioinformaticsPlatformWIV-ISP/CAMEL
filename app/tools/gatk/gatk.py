import logging
import os
import re
import abc

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class GATK(Tool):

    """
    Super class for GATK tools
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, tool_name, version, camel):
        """
        Initialize a GATK tool
        :param tool_name: Tool name
        :param version: Tool version
        :param camel: Camel instance
        :return: None
        """
        super(GATK, self).__init__(tool_name, version, camel)
        self._function_name = ''
        self._specific_parameters = []
        self._required_inputs = []
        self._input_string = ''
        self._output_string = ''
        self._option_string = ''
        self._output_type = ''

    def _execute_tool(self):
        """
        Run a GATK function
        :return: None
        """
        self._set_input()
        self._set_output()
        self._set_specific_parameters()
        self._build_command()
        self._execute_command()
        self._set_informs()

    def _check_input(self):
        """
        Check input for a tool and prepare command line parameters for input
        :return: None
        """
        for input_file in self._required_inputs:
            if input_file not in self._tool_inputs:
                raise InvalidInputSpecificationError('GATK {!r} required {!r} input is missing in _tool_inputs!'.format(
                    self._function_name, input_file))

    def _set_input(self):
        """
        Set the input specification
        :return: None
        """
        if 'TXT_intervals' in self._tool_inputs:
            self._input_string += "-L {} ".format(self._tool_inputs['TXT_intervals'][0].path)

        if 'VCF' in self._tool_inputs:
            self._input_string += "-V {} ".format(self._tool_inputs['VCF'][0].path)

        if 'FASTA_REF' in self._tool_inputs:
            self._input_string += "-R {} ".format(self._tool_inputs['FASTA_REF'][0].path)

    def _set_output(self):
        """
        Set the output specification
        :return: None
        """
        self._tool_outputs[self._output_type] = [
            ToolIOFile(os.path.join(self._folder, self._parameters['output'].value))
        ]

    def _set_specific_parameters(self):
        """
        Set specific parameters that need special handling when required
        :return: None
        """
        pass

    def _build_command(self):
        """
        Build the command to run tool
        :return: None
        """
        self._option_string += " ".join(self._build_options(excluded_parameters=self._specific_parameters))
        self._command.command = " ".join([
            self._tool_command, self._input_string, self._output_string, self._option_string
        ])

    def _check_command_output(self):
        """
        Check the result of GATK tool run
        :return: None
        """
        if not re.match('Exit status: 0', self.stdout.split('\n')[-2].rstrip()):
            raise ToolExecutionError("GATK tool {} fails to run, message: \n{}".format(
                self._function_name, self.stdout))

    def _set_informs(self):
        """
        Set informs by analyze the output
        :return: None
        """
        # log WARNINGS in info.log
        for l in self.stdout.split('\n'):
            if re.match('WARNING', l):
                logging.info(" GATK - {}".format(l))

            # E.g., The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428
            match = re.search('The Genome Analysis Toolkit (GATK) (.+),', l)
            if match:
                self.informs['tool_name'] = 'GATK'
                self.informs['tool_version'] = match.group(1)
            # Total filtering statistics
            #   0 reads were filtered out during the traversal out of approximately 1090654 total reads (0.00%)
            match = re.search(
                '(\d+) reads were filtered out during the traversal out of approximately (\d+) total reads \((.+%)\)', l)
            if match:
                self.informs['reads_total'] = match.group(2)
                if match.group(1) != '0':
                    self.informs['filtered_reads'] = "{}({})".format(match.group(1), match.group(3))
            # per Filter statistics:
            #    0 reads (0.00% of total) failing BadCigarFilter
            match = re.search('(\d+) reads \((.+%) of total\) failing (.+)', l)
            if match and match.group(1) != '0':
                self.informs[match.group(3)] = "{}({})".format(match.group(1), match.group(2))
