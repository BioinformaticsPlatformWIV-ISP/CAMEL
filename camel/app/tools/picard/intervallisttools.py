import glob

from camel.app.camel import Camel
from camel.app.tools.picard.picard import Picard
from camel.app.io.tooliofile import ToolIOFile

class IntervalListTools(Picard):
    """
    Class for Picard IntervalListTools function
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard IntervalListTools', '2.23.3', camel)

        self._supported_inputs = ['TXT_intervalListFiles']

    def _check_input(self) -> None:
        """
        Set the input specification, this default function handles only one SAM or BAM file as input
        :return: None
        """
        super(Picard, self)._check_input()

        self._set_input()

    def _set_input(self) -> None:
        """
        Set the input specification
        Overrides method in parent class.
        :return: None
        """
        if 'TXT_intervals' in self._tool_inputs:
            self._input_string += f" INPUT={self._tool_inputs['TXT_intervals'][0].path}"

    def _set_output(self) -> None:
        """
        Set the output specification
        Overrides method in parent class.
        :return: None
        """
        interval_lists = glob.glob(self._parameters['output_dir'].value + "/temp_*/scattered.interval_list")
        self._tool_outputs['TXT_intervalLists'] = []

        for list in interval_lists:
            self._tool_outputs['TXT_intervalLists'].append(ToolIOFile(list))

    def _execute_tool(self) -> None:
        """
        Function to run Picard function.
        Overrides method in parent class.
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._set_output()
        self._set_informs()
