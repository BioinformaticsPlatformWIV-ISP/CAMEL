import os

from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.tools.velvet.velvet import Velvet


class Velveth(Velvet):
    """
    Velveth generates kmer hash database from reads data to be assembled by velvetg
    """

    def __init__(self, camel):
        """
        Initialize Velveth
        :param camel: Camel instance
        :return: None
        """
        super().__init__('velveth', '1.2.10', camel)

    def _execute_tool(self):
        """
        Function to run BWA index
        :return: None
        """
        self._set_input()
        self.__build_command()
        self._execute_command()
        self.__set_output()
        self.__set_inform()

    def __build_command(self):
        """
        Build the command to run tool
        :return: None
        """
        self._command.command = "{} {} {} {} {}".format(
            self._tool_command,
            self._parameters['output_dir'].value,
            self._parameters['kmer'].value,
            " ".join(self._build_options(excluded_parameters=['kmer', 'output_dir'])),
            self._input_string
        )

    def __set_output(self):
        """
        Specify the output of tool and the command line options
        :return: None
        """
        self._tool_outputs['DIR_DB'] = [
            ToolIODirectory(os.path.join(self._folder, self._parameters['output_dir'].value))
        ]

    def __set_inform(self):
        """
        Set the tool inform
        :return: None
        """
        self.informs['kmer'] = self._parameters['kmer'].value
