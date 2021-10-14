from pathlib import Path

from camel.app.camel import Camel
from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4MarkDuplicatesSpark(GATK4):

    """
    Class for GATK MarkDuplicatesSpark function
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize a GATK MarkDuplicatesSpark tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('gatk4 MarkDuplicatesSpark', '4.1.9.0', camel)

        self._required_inputs = ['BAM']
        self._output_type = 'BAM'
        self._specific_parameters = ['threads']

    def _set_input(self) -> None:
        """
        Set the input specification
        :return: None
        """
        super(GATK4MarkDuplicatesSpark, self)._set_input()

        bam_file = self._tool_inputs['BAM'][0].path
        self._input_string += f"-I {bam_file} "

    def _set_specific_parameters(self) -> None:
        """
        Sets the parameters for the number of threads and the temporary directory
        to use. The temporary directory is ideally a locally attached SSD disk.
        :return: None
        """
        self._option_string = f'--conf \'spark.executor.cores={self._parameters["threads"].value}\' '
        if Path('/scratch').is_dir():
            self._option_string += f'--conf \'spark.local.dir=/scratch\' '
        elif Path('/temp').is_dir():
            self._option_string += f'--conf \'spark.local.dir=/temp\' '
