from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.tools.tool import Tool


class characterize_neisseria_capsule(Tool):
    """
    characterize_neisseria_capsule is a tool implementing a WGS-based method for N. meningitidis
    serogroup prediction. it identifies capsule genes and genetic variations that might impact their expression.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the characterize_neisseria_capsule tool.
        :param camel: CAMEL instance
        """
        super().__init__('characterize_neisseria_capsule', 'a75a009', camel)

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - FASTA_dir is the only required input
        :return: None
        """
        if 'FASTA_dir' not in self._tool_inputs:
            raise InvalidInputSpecificationError(
                f'Required input directory containing the fasta files is missing')
        if 'OUTPUT_dir' not in self._tool_inputs:
            raise InvalidInputSpecificationError(
                f'Required output directory path is missing')
        if 'THREADS' not in self._tool_inputs:
            raise InvalidInputSpecificationError(
                f'Number of threads to use is missing')
        super()._check_input()

    def __build_input_string(self) -> str:
        """
        Creates the string with the input files
        :return: String with the input parameters
        """
        inputs = []
        if 'FASTA_dir' in self._tool_inputs:
            inputs.append(f"-d {self._tool_inputs['FASTA_dir'][0].path}")
        if 'OUTPUT_dir' in self._tool_inputs:
            inputs.append(f"-o {self._tool_inputs['OUTPUT_dir'][0].path}")
        if 'THREADS' in self._tool_inputs:
            inputs.append(f"-t {self._tool_inputs['THREADS'][0].path}")
        return ' '.join(inputs)

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command.
        :return: None
        """
        input_string = self.__build_input_string()
        self._command.command = ' '.join([self._tool_command, input_string])

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
