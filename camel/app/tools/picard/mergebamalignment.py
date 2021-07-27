import re

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.tools.picard.picard import Picard


class MergeBamAlignment(Picard):

    """
    Class for Picard MergeBamAlignment function
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard MergeBamAlignment', '2.23.3', camel)

        self._function_name = 'MergeBamAlignment'
        self._main_inputs = []
        self._extra_inputs = ['BAM_UNMAPPED', 'BAM_ALIGNED', 'FASTA_REF']
        self._specific_parameters = ["attributes_to_remove_multi"]

    def _check_input(self) -> None:
        """
        Check for required input files
        :return: None
        """
        if 'BAM_UNMAPPED' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Picard MergeBamAlignment: Input file BAM_UNMAPPED is not defined")

        if 'FASTA_REF' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Picard MergeBamAlignment: Input file FASTA_REF is not defined")

    def _set_input(self) -> None:
        """
        Set required inputs specification
        :return: None
        """
        self._input_string += f"UNMAPPED={self._tool_inputs['BAM_UNMAPPED'][0].path} "

        self._input_string += f"REFERENCE_SEQUENCE={self._tool_inputs['FASTA_REF'][0].path} "

        if 'BAM_ALIGNED' in self._tool_inputs:
            self._input_string +=  f"ALIGNED={self._tool_inputs['BAM_ALIGNED'][0].path} "

    def _build_command(self) -> None:
        """
        Build the command to run tool
        :return: None
        """
        build_options = self._build_options(excluded_parameters=self._specific_parameters, delimiter='=')

        if 'attributes_to_remove_multi' in self._parameters:
            attributes = self._parameters['attributes_to_remove_multi'].value.split(",")
            for attribute in attributes:
                build_options.append(f"ATTRIBUTES_TO_REMOVE={attribute}")

        option_string = " ".join(build_options)

        self._command.command = " ".join([
            "java", self._java_options, "-jar $PICARD_JAR", self._tool_command, self._java_options_temp_dir, self._input_string, self._output_string,
            option_string, '2>&1'
        ])

    def _set_informs(self) -> None:
        """
        Analyse the result of picard run and update tool.informs
        :return: None
        """
        for line in self.stdout.splitlines():
            m = re.search(r'Finished reading (\d+) total records from alignment SAM/BAM.', line)
            if m:
                self.informs['reads_total'] = m.group(1)
            m = re.search(r'Wrote (\d+) alignment records and (\d+) unmapped reads.', line)
            if m:
                if m.group(1):
                    self.informs['aligned_reads'] = m.group(1)
                if m.group(2):
                    self.informs['unmapped_reads'] = m.group(2)
