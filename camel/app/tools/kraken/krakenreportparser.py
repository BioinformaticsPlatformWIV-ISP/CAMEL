import logging

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.tools.tool import Tool


class KrakenReportParser(Tool):
    """
    Parses Kraken output reports.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Kraken Report Parser', '0.1', camel)

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self._informs['contaminants_warn'] = []
        self._informs['contaminants_fail'] = []
        with open(self._tool_inputs['TSV'][0].path) as handle:
            for line in handle.readlines():
                parts = line.split('\t')
                if parts[3] != 'S':
                    continue
                percentage = float(parts[0])
                if percentage < float(self._parameters['threshold_warn'].value):
                    continue
                species_name = parts[-1].strip()
                if species_name == self._parameters['expected_species'].value:
                    self._informs['expected'] = (species_name, percentage)
                elif percentage < float(self._parameters['threshold_fail'].value):
                    self._informs['contaminants_warn'].append((species_name, percentage,))
                else:
                    self._informs['contaminants_fail'].append((species_name, percentage,))

        if 'expected' not in self._informs:
            logging.warning("No reads matching the expected species found!")
            self._informs['expected'] = (self._parameters['expected_species'].value, 0)

        self._informs['contaminants_warn'].sort(key=lambda x: -x[1])
        self._informs['contaminants_fail'].sort(key=lambda x: -x[1])

    def _check_input(self):
        """
        Checks the tool input.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("TSV input is required.")
        super()._check_input()
