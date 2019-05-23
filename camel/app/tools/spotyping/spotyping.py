import json
import logging

import os

from camel.app.command.command import Command
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class SpoTyping(Tool):
    """
    SpoTyping: fast and accurate in silico Mycobacterium spoligotyping from sequence reads.

    Input:
        - FASTQ: 1 (SE) or 2 (PE) FASTQ files

    Output:
        - VAL_type_binary: Binary spoligotype
        - VAL_type_octal: Octal spoligotype
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('SpoTyping', '2.1', camel)

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTQ' not in self._tool_inputs:
            raise InvalidInputSpecificationError("FASTQ input is required")
        elif not (0 < len(self._tool_inputs['FASTQ']) <= 2):
            raise InvalidInputSpecificationError("Only 1 (SE) or 2 (PE) FASTQ inputs are supported")
        super()._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join([f.path for f in self._tool_inputs['FASTQ']]),
            ' '.join(self._build_options())
        ])
        self._execute_command()
        type_binary, type_octal = self._parse_output_file()
        self._tool_outputs['VAL_type_binary'] = [ToolIOValue(type_binary)]
        self._tool_outputs['VAL_type_octal'] = [ToolIOValue(type_octal)]
        self._tool_outputs['LOG'] = [ToolIOFile(os.path.join(self._folder, '{}.log'.format(
            self._parameters['output_basename'].value)))]

    def _parse_output_file(self):
        """
        Parses the output file.
        :return: Spoligotype (Binary), Spoligotype (Octal)
        """
        output_file_path = os.path.join(self._folder, self._parameters['output_basename'].value)
        if not os.path.isfile(output_file_path):
            raise ToolExecutionError("Output file not found")
        with open(output_file_path, 'r') as handle:
            try:
                _, type_binary, type_octal = handle.readlines()[-1].strip().split('\t')
                return type_binary, type_octal
            except IndexError:
                raise ToolExecutionError("Output file has an invalid format")

    def _check_command_output(self):
        """
        Checks the command output to checks if the tool executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            last_line = self._command.stderr.splitlines()[-1]
            if last_line.startswith('urllib2.URLError'):
                logging.warning('Could not contact SITVIT server')
            else:
                raise ToolExecutionError(last_line)

    def _extract_metadata(self, type_octal):
        """
        Extracts the metadata for the detected Spoligotype.
        :return: Spoligotype metadata
        """
        command = Command('{} echo $SPOTYPING_METADATA'.format(self._build_dependencies()))
        command.run_command(self._folder)
        metadata_path = command.stdout.strip()
        with open(metadata_path) as handle:
            metadata = json.load(handle)
        keys = ('SIT', 'geo', 'label', 'total')
        if type_octal in metadata:
            return {k: metadata[type_octal][k] for k in keys}
        else:
            return {k: 'NA' for k in keys}
