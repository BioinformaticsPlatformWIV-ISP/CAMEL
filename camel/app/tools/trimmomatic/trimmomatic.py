import re

import os

from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Trimmomatic(Tool):

    """
    A flexible read trimming tool for Illumina NGS Data.
    """

    def __init__(self, camel):
        """
        Initializes Trimmomatic.
        :param camel: Camel instance
        """
        super().__init__('Trimmomatic', '0.38', camel)
        self._mode = None

    def _execute_tool(self):
        """
        Runs Trimmomatic.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()
        self.__set_informs()
        if not self._informs.get('succeed', 'False'):
            raise ToolExecutionError("Error running trimmomatic")

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if 'FASTQ_PE' in self._tool_inputs:
            self._mode = 'PE'
            if len(self._tool_inputs['FASTQ_PE']) != 2:
                raise ValueError("Paired end input requires exactly 2 files.")
        elif 'FASTQ_SE' in self._tool_inputs:
            self._mode = 'SE'
            if len(self._tool_inputs['FASTQ_SE']) != 1:
                raise ValueError("Single end input requires exactly 1 file.")
        else:
            raise ValueError("No FASTQ_PE of FASTQ_SE input found")
        super(Trimmomatic, self)._check_input()

    def __build_command(self):
        """
        Builds the command.
        :return: None
        """
        if self._mode == 'PE':
            options = self.__build_pe_command()
        else:
            options = self.__build_se_command()
        options += self._build_options(excluded_parameters=['baseout', 'threads', 'illuminaclip_PE', 'illuminaclip_SE'],
                                       delimiter='')
        self._command.command = '{} {}'.format(self._tool_command, ' '.join(options))

    def __build_se_command(self):
        """
        Builds the command to run in single end mode.
        :return: Command options
        """
        options = [self._mode]
        if 'threads' in self._parameters:
            options.append(str(self._parameters['threads']))
        options.append(self._tool_inputs['FASTQ_SE'][0].path)
        if 'baseout' in self._parameters:
            options.append(str(self._parameters['baseout']))
        if 'illuminaclip_SE' in self._parameters:
            options.append(self._parameters['illuminaclip_SE'].option + self._parameters['illuminaclip_SE'].value)

        return options

    def __build_pe_command(self):
        """
        Builds the command to run in paired end mode.
        :return: Command options
        """
        options = [self._mode]
        if 'baseout' in self._parameters:
            options.append(str(self._parameters['baseout']))
        if 'threads' in self._parameters:
            options.append(str(self._parameters['threads']))
        options.append(' '.join(f.path for f in self._tool_inputs['FASTQ_PE']))
        if 'illuminaclip_PE' in self._parameters:
            options.append(self._parameters['illuminaclip_PE'].option + self._parameters['illuminaclip_PE'].value)

        return options

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        basename = os.path.splitext(self._parameters['baseout'].value)[0]
        if self._mode == 'PE':
            self._tool_outputs['FASTQ_PE'] = [
                ToolIOFile(os.path.join(self._folder, basename + '_1P.fastq')),
                ToolIOFile(os.path.join(self._folder, basename + '_2P.fastq'))
            ]
            self._tool_outputs['FASTQ_SE_FORWARD'] = [ToolIOFile(os.path.join(self._folder, basename + '_1U.fastq'))]
            self._tool_outputs['FASTQ_SE_REVERSE'] = [ToolIOFile(os.path.join(self._folder, basename + '_2U.fastq'))]
        else:
            self._tool_outputs['FASTQ'] = [ToolIOFile(os.path.join(self._folder, basename + '.fastq'))]
        self.__remove_empty_outputs()

    def __remove_empty_outputs(self):
        """
        Removes the empty files from the outputs.
        :return: None
        """
        for key in self._tool_outputs:
            for tool_output_file in self._tool_outputs[key]:
                if tool_output_file.size == 0:
                    self._tool_outputs[key].remove(tool_output_file)

    def __set_informs(self):
        """
        Adds the trimming statistics to the informs.
        :return: None
        """
        self._informs['mode'] = self._mode
        for line in self.stdout.splitlines():
            qc_encoding = re.search("Quality encoding detected as (?P<encode>\\w+)", line)
            if qc_encoding:
                self._informs['encoding'] = qc_encoding.group('encode')
            elif re.match("Input Read Pairs", line):
                res = re.match(r"Input Read Pairs: (\d+) Both Surviving: (\d+ \([\d.%]+\)) "
                               r"Forward Only Surviving: (\d+ \([\d.%]+\)) "
                               r"Reverse Only Surviving: (\d+ \([\d.%]+\)) "
                               r"Dropped: (\d+ \([\d.%]+\))", line)
                self._informs['paired_reads_in'] = res.groups()[0]
                self._informs['paired_reads_out'] = res.groups()[1]
                self._informs['forward_only_reads'] = res.groups()[2]
                self._informs['reverse_only_reads'] = res.groups()[3]
                self._informs['reads_drop'] = res.groups()[4]
            elif re.match('Input Read', line):
                res = re.match(r"Input Reads: (\d+) "
                               r"Surviving: (\d+ \([\d.%]+\)) "
                               r"Dropped: (\d+ \([\d.%]+\))", line)
                self._informs['reads_in'] = res.groups()[0]
                self._informs['reads_out'] = res.groups()[1]
                self._informs['reads_drop'] = res.groups()[2]
            elif line.strip() == 'Exit status: 0':
                self._informs['succeed'] = True
