import os
import re

from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class FastQC(Tool):

    """
    FastQC tool.
    """

    def __init__(self, camel):
        """
        Initializes FastQC.
        :param camel: Camel instance
        """
        super(FastQC, self).__init__('FastQC', '0.11.5', camel)

    def _execute_tool(self):
        """
        Runs FastQC.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if 'FASTQ' not in self._tool_inputs or len(self._tool_inputs['FASTQ']) == 0:
            raise ValueError("Required FASTQ input file is missing for FastQC.")
        super(FastQC, self)._check_input()

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command to run fastQC
        :return: none
        """
        self._command.command = ' '.join([self._tool_command,
                                          ' '.join(in_file.path for in_file in self._tool_inputs['FASTQ']),
                                          '--outdir .',
                                          ' '.join(self._build_options())])

    def _check_command_output(self):
        """
        Checks if the command output is valid.
        :return: None
        """
        for line in self.stderr.splitlines():
            if not line.startswith('Started analysis'):
                raise ToolExecutionError("Error executing FastQC: {}".format(self.stderr.strip()))

    @staticmethod
    def __get_output_folder(execution_folder, input_file):
        """
        Returns the output folder for the given input file.
        :param execution_folder: Folder where the command is executed
        :param input_file: Input file name
        :return: Output folder
        """
        sample_base_name = input_file.basename.split('.')[0]
        if re.search(r'\.fastq$', sample_base_name):
            sample_base_name = re.sub(r'\.fastq$', '', sample_base_name)
        for sub_folder in os.listdir(execution_folder):
            if sub_folder.startswith(sample_base_name) and sub_folder.endswith('_fastqc'):
                full_path = os.path.join(execution_folder, sub_folder)
                if os.path.isdir(full_path):
                    return full_path
        raise IOError("No output directory for FastQC input {} found.".format(input_file))

    @staticmethod
    def _analyze_summary_file(summary_file):
        """
        Analyze fastqc output summary.txt (of a given input file)
        :param summary_file: FastQC summary file
        :return: Dictionary containing the summary information
        """
        summary_info = {'passed': True, 'warnings': [], 'fails': []}
        with open(summary_file, 'r') as input_handle:
            for line in input_handle.readlines():
                status, test_name, _ = line.split('\t')
                if status == 'WARN':
                    summary_info['warnings'].append(test_name)
                elif status == 'FAIL':
                    summary_info['fails'].append(test_name)
        if len(summary_info['fails']) > 0:
            summary_info['passed'] = False
        return summary_info

    def __set_output(self):
        """
        Set the output of FastQC.
        :return: None
        """
        self._tool_outputs['HTML'] = []
        self._tool_outputs['TXT'] = []
        for input_file in self._tool_inputs['FASTQ']:
            output_folder = self.__get_output_folder(self._folder, input_file)
            self.informs[input_file.path] = FastQC._analyze_summary_file(os.path.join(output_folder, 'summary.txt'))
            for output_file in os.listdir(output_folder):
                if output_file == 'fastqc_report.html':
                    self._tool_outputs['HTML'].append(ToolIOFile(os.path.join(output_folder, output_file)))
                elif output_file == 'fastqc_data.txt':
                    self._tool_outputs['TXT'].append(ToolIOFile(os.path.join(output_folder, output_file)))
