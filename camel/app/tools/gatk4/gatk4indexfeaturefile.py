from camel.app.camel import Camel
from camel.app.tools.gatk4.gatk4 import GATK4
from camel.app.io.tooliofile import ToolIOFile
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError


class GATK4IndexFeatureFile(GATK4):
    """
    ==============================
    GATK IndexFeatureFile 4.1.9.0
    ==============================
    Creates an index file for the various kinds of feature-containing files supported by GATK

    Required inputs:
    ----------------
    'VCF'|'BED':        ToolIOFile object. Input VCF or BED file.

    Output (optional):
    -------
    'IDX':              ToolIOFile object. Output index file.  If missing, the tool will create an
                        index file in the same directory as the input file.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize GATK4IndexFeatureFile tool.
        :param camel: Camel instance
        :return: None
        """
        super(GATK4IndexFeatureFile, self).__init__('gatk4 GATK4IndexFeatureFile', '4.1.9.0', camel)

        self._required_inputs = ['VCF', 'VCF_gz', 'BED']
        self._input_format = []
        self._output_type = 'IDX'

    def _check_input(self) -> None:
        """
        Check input for the tool. Checks whether at least on of the required inputs is specified
        :return: None
        """
        for input_file in self._required_inputs:
            if input_file in self._tool_inputs:
                self._input_format.append(input_file)

        if len(self._input_format) != 1:
            raise InvalidInputSpecificationError(''.join(f'GATK {self._name} requires either one of: ', ' '.join(self._required_inputs)))

        super(GATK4, self)._check_input()

    def _set_input(self) -> None:
        """
        Set input for a tool
        :return: None
        """
        self._input_string += f"-I {self._tool_inputs[self._input_format[0]][0].path} "

    def _build_command(self) -> None:
        """
        Build the command to run tool
        :return: None
        """
        super(GATK4IndexFeatureFile, self)._build_command()

        if 'output' in self._parameters:
            self._command.command += f"-o {self._parameters['output'].value}"

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        if 'output' in self._parameters:
            self._tool_outputs[self._output_type] = [
                ToolIOFile(self.folder / self._parameters['output'].value)
            ]
        else:
            if self._input_format[0] == 'VCF_gz':
                self._tool_outputs[self._output_type] = [
                    ToolIOFile(self.folder / f"{self._tool_inputs[self._input_format[0]][0].path}.tbi")]
            else:
                self._tool_outputs[self._output_type] = [
                    ToolIOFile(self.folder / f"{self._tool_inputs[self._input_format[0]][0].path}.idx")]
