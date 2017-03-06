import os

from app.io.tooliofile import ToolIOFile
from app.tools.samtools.samtools import Samtools


class SamtoolsDepth(Samtools):
    """
    Calculates the coverage depth of an alignment.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super(SamtoolsDepth, self).__init__('samtools depth', '1.3.1', camel)

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if 'BAM' not in self._tool_inputs:
            raise ValueError("No BAM input file found")
        if len(self._tool_inputs['BAM']) != 1:
            raise ValueError("Exactly one BAM input file expected")
        super(Samtools, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __build_command(self):
        """
        Builds the command.
        :return: None
        """
        self._command.command = ' '.join(
            [self._tool_command,
             ' '.join(self._build_options(['output_filename'])),
             self._tool_inputs['BAM'][0].path,
             ' > {}'.format(self._parameters['output_filename'].value)])

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        output_file_path = os.path.join(self._folder, self._parameters['output_filename'].value)
        self._tool_outputs['TXT'] = [ToolIOFile(output_file_path)]
        self._informs['median_depth'] = SamtoolsDepth.calculate_median_coverage(output_file_path)

    @staticmethod
    def median(input_list):
        """
        Returns the median value of a list.
        :return:
        """
        sorted_list = sorted(input_list)
        middle = len(input_list) // 2
        if len(input_list) % 2:
            return sorted_list[middle]
        else:
            median = (sorted_list[middle] + sorted_list[middle - 1]) / 2
            return median

    @staticmethod
    def calculate_median_coverage(output_path):
        """
        Calculates the median coverage.
        :param output_path: Path to the output files.
        :return: None
        """
        coverage_values = []
        with open(output_path) as output_file:
            for line in output_file.readlines():
                seq_id, pos, count = line.split('\t')
                coverage_values.append(int(count))
        return SamtoolsDepth.median(coverage_values)
