import abc

from app.tools.tool import Tool


class Seqtk(Tool):
    """
    Base class for all seqtk functionality
    """

    def __init__(self, tool_name, version, camel):
        """
        Initialize seqtk
        :param tool_name: Tool name
        :param version: Tool version
        :param camel: Camel instance
        :return: None
        """
        super(Seqtk, self).__init__(tool_name, version, camel)
        self._function_name = ''
        # parameters that should not be handled by self.build_options function
        self._specific_parameters = []
        # alternative types of files that can be used as main input of a seqtk tool
        # - FASTA or FASTQ
        self._supported_inputs = ['FASTA', 'FASTQ', 'FASTA_SE', 'FASTQ_SE', 'FASTA_PE', 'FASTQ_PE']
        self.input_mode = None       # 'SE' or 'PE'
        self.input_file_type = None  # 'FASTQ' or 'FASTA'
        self._input_string = ''
        self._output_string = ''
        self._input_files = []
        self._output_files = []

    def _execute_tool(self):
        """
        Function to run BWA index
        :return: None
        """
        self._set_output()
        self._build_command()
        self._execute_command()

    def _check_input(self):
        """
        Check input requirements to run a seqtk tool
        :return: None
        """
        if len(self._supported_inputs) != 0:
            input_type = self.__get_supported_input_type()
            self._input_files = [f.path for f in self._tool_inputs[input_type]]
            self.__check_supported_input_files()

    def __get_supported_input_type(self):
        """
        Check the type of supported_inputs (alternatives but still required) specified in _tool_inputs and return it
        :return: type of supported_input found
        """
        for input_type in self._supported_inputs:
            if input_type in self._tool_inputs:
                type_inform = input_type.split("_")
                if len(type_inform) > 1:
                    self.input_mode = type_inform[1]
                else:
                    # by default, support only one input
                    self.input_mode = 'SE'
                self.input_file_type = type_inform[0]
                return input_type

        raise KeyError('Seqtk function {!r} required {!r} input is missing!'.format(
            self._function_name, self._supported_inputs))

    def __check_supported_input_files(self):
        """
        Check supported input files are correct
        :return: None
        """
        if self.input_mode == 'SE' and len(self._input_files) != 1:
            raise ValueError("Seqtk function {} SE mode supports only one input file.".format(self._function_name))
        elif self.input_mode == 'PE' and len(self._input_files) != 2:
            raise ValueError("Seqtk function {} PE mode supports only two input files.".format(self._function_name))

    @abc.abstractmethod
    def _set_output(self):
        """
        Set the output specification
        :return: None
        """
        return

    def _build_command(self):
        """
        Build the command to run tool
        :return: None
        """
        self._command.command = "{} {} {} > {}".format(
            self._tool_command,
            " ".join(self._build_options(excluded_parameters=self._specific_parameters)),
            self._input_string,
            self._output_string
        )
