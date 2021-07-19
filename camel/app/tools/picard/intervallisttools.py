from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.picard.picard import Picard


class IntervalListTools(Picard):
    """
    ==============================
    Picard IntervalListTools 2.23.3
    ==============================
    This tool offers multiple interval list file manipulation capabilities, including: sorting, merging, subtracting,
    padding, and other set-theoretic operations.

    Required inputs:
    ----------------
    'TXT_intervalListFiles' ToolIOFile object. One or more interval lists. Supported formats are interval_list and VCF.
                            IntervalList should be denoted with the extension .interval_list, while a VCF must have one
                            of .vcf, .vcf.gz, .bcf

    Output:
    -------
    'TXT_intervalListFiles' ToolIOFile object. One or more interval lists

    """
    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard IntervalListTools', '2.23.3', camel)

        self._supported_inputs = ['TXT_intervalListFiles', 'VCF']

    def _check_input(self) -> None:
        """
        Check the input specification
        Overrides method in parent class
        :return: None
        """
        # Method from parent class of the parent class
        super(Picard, self)._check_input()

        self._set_input()

    def _set_input(self) -> None:
        """
        Set the input specification
        Overrides method in parent class.
        :return: None
        """
        self._input_string = ""

        # one or more input files possible
        for io_interval_file in self._tool_inputs['TXT_intervals']:
            self._input_string += f" INPUT={io_interval_file.path}"

    def _set_output(self) -> None:
        """
        Set the output specification
        Overrides method in parent class.
        :return: None
        """
        scatter_count = int(self._parameters.get('scatter_count', 1).value)

        if scatter_count == 1:
            self._tool_outputs['TXT_intervalLists'] = [ToolIOFile(Path(self._folder) / self._parameters['output'].value)]
        else:
            intervallists_list = []
            for n in range(1, scatter_count + 1):
                intervallists_list.append(ToolIOFile(Path(self._folder) / self._parameters['output'].value / 'temp_{:04d}_of_{:d}'.format(n, scatter_count) / "scattered.interval_list" ))
            self._tool_outputs['TXT_intervalLists'] = intervallists_list

