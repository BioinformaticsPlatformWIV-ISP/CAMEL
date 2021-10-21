from pathlib import Path

from camel.app.camel import Camel
from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4GenomicsDBImport(GATK4):
    """
    =============================
    GATK GenomicsDBImport 4.1.9.0
    =============================
    Import single-sample GVCFs into GenomicsDB before joint genotyping.

    Required inputs:
    ----------------
    'gVCF':             ToolIOFile object. Multiple gVCF files.
    'TXT_intervals':    ToolIOFile object. Contains a gatk interval list (.list or .intervals)

    Output:
    -------
    None

    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize the GenomicsDBImport tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('gatk4 GenomicsDBImport', '4.1.9.0', camel)

        self._required_inputs = ['gVCF', 'TXT_intervals']
        self._output_type = 'TXT_success'
        self._specific_parameters = ['output']
        self._java_options = f'"-mx8G -ms8G -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Djava.io.tmpdir={self._temp_dir}"'

    def _set_input(self) -> None:
        """
        Set the input specification
        Overrides method in parent class.
        :return: None
        """
        super(GATK4GenomicsDBImport, self)._set_input()

        for f in self._tool_inputs['gVCF']:
            self._input_string += f"--variant {f.path} "

    def _check_command_output(self) -> None:
        """
        Check the result of the GATK run
        :return: None
        """
        super(GATK4GenomicsDBImport, self)._check_command_output()

        # Snakemake placeholder
        Path.touch(Path(self._parameters['output'].value))