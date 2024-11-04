from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.star.star import Star


class StarAlign(Star):
    """
    align spliced transcripts with STAR
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes STAR 2.7.11b alignment
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('STAR', '2.7.11b', camel)
        self._required_inputs = ['FASTQ', 'INDEX_DIR']
        self._input_string = "--runMode alignReads"

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if len(self._tool_inputs['FASTQ']) > 2:
            raise ValueError("Too many inputs of type 'FASTQ'")
        super()._check_input()

    def _set_input(self) -> None:
        """
        Set the input specification and the input string
        :return: None
        """
        self._input_string += " --readFilesIn"
        for fastq in self._tool_inputs['FASTQ']:
            self._input_string += f" {Path(str(fastq))}"

        if 'GTF' in self._tool_inputs:
            self._input_string += f" --sjdbGGTFfile {Path(str(self._tool_inputs['GTF'][0]))}"

        self._input_string += f" --genomeDir {Path(str(self._tool_inputs['INDEX_DIR'][0]))}"

    def _set_output(self) -> None:
        """
        Set the output specification and the output string
        :return: None
        """
        filename_output = self._parameters['filename_output'].value
        self._output_string += f" > {filename_output}"
        self._tool_outputs['ALIGNMENT'] = [ToolIOFile(self.folder / filename_output)]
