import json
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Cutadapt(Tool):
    """
    =============================
    Cutadapt 3.5
    =============================
    Cutadapt trims adapters from high-throughput sequencing reads

    Input:
    ------
    FASTQ_PE|FASTQ_SE:      ToolIOFile object. Two resp. one FASTQ files.

    Output:
    -------
    FASTQ_PE|FASTQ_SE:      ToolIOFile object. Two resp. one FASTQ trimmed files.
    JSON_report:            ToolIOFile object. JSON file containing all stats of the run. The contents of this file
                            are added to the tool informs

    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('cutadapt', '3.5', camel)
        self._input_string = ''
        self._output_string = ''

    def _execute_tool(self) -> None:
        """
        Run a cutadapt function
        :return: None
        """
        self._set_input()
        self._set_output()
        self._build_command()
        self._execute_command()
        self._set_informs()

    def _check_input(self) -> None:
        """
        Check input of cutadapt function. Overrides method in parent class
        :return: None
        """
        super(Cutadapt, self)._check_input()

        if 'FASTQ_PE' in self._tool_inputs:
            if len(self._tool_inputs['FASTQ_PE']) != 2:
                raise ValueError("Paired end input requires exactly 2 files.")
        elif 'FASTQ_SE' not in self._tool_inputs:
            raise ValueError("No FASTQ_PE or FASTQ_SE input found")

    def _set_input(self) -> None:
        """
        Set input of cutadapt function
        :return: None
        """
        if 'FASTQ_PE' in self._tool_inputs:
            self._input_string = " ".join([str(f.path) for f in self._tool_inputs['FASTQ_PE']])
        elif 'FASTQ_SE' in self._tool_inputs:
            self._input_string = str(self._tool_inputs['FASTQ_SE'][0].path)

    def _set_output(self) -> None:
        """
        Set output of cutadapt function
        :return: None
        """
        # Set output string for command
        if 'FASTQ_PE' in self._tool_inputs:
            self._output_string = f"-o {Path(self.folder) / self._parameters['output_basename'].value}_R1_trimmed.fastq.gz " \
                            f"-p {Path(self.folder) / self._parameters['output_basename'].value}_R2_trimmed.fastq.gz "
        elif 'FASTQ_SE' in self._tool_inputs:
            self._output_string = f"-o {Path(self.folder) / self._parameters['output_basename'].value}_trimmed.fastq.gz "

        # Set tool output
        if 'FASTQ_PE' in self._tool_inputs:
            self._tool_outputs['FASTQ'] = [ToolIOFile(Path(self.folder) / f"{self._parameters['output_basename'].value}_R1_trimmed.fastq.gz"),
                                           ToolIOFile(Path(self.folder) / f"{self._parameters['output_basename'].value}_R2_trimmed.fastq.gz")]
        elif 'FASTQ_SE' in self._tool_inputs:
            self._tool_outputs['FASTQ'] = [ToolIOFile(Path(self.folder) / f"{self._parameters['output_basename'].value}_trimmed.fastq.gz")]

        if self._parameters['report_json']:
            self._tool_outputs['JSON_report'] = [ToolIOFile(Path(self.folder) / self._parameters['report_json'].value)]

    def _build_command(self) -> None:
        """
        Build command of cutadapt function
        :return: None
        """
        self._command.command = ' '.join([self._tool_command,
                                          *self._build_options(excluded_parameters=['output_basename']),
                                          self._output_string,
                                          self._input_string])

    def _set_informs(self) -> None:
        """
        Set informs for cutadapt function
        :return: None
        """
        with open(self.folder / self._parameters["report_json"].value) as json_file:
            data = json.load(json_file)
            self._informs.update(data)
        