from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.loggers import logger
from camel.app.tools.tool import Tool

class KrakenReportParserFasta(Tool):
    """
    Parses Kraken output reports for fasta input to the species (S) level.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Kraken Report Parser Fasta', '0.1', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        allowed_species = (
            self._parameters['allowed_species'].value.split(',')) if 'allowed_species' in self._parameters else []
        if self._parameters['level_of_depth'].value not in ['S', 'G']:
            logger.error('Please choose either "G" or "S" for parameter level_of_depth')

        self._informs['contaminants_fail'] = []
        self._informs['allowed'] = []
        self._informs['contaminants_warn'] = []

        with open(self._tool_inputs['TSV_out'][0].path) as report:
            total_bp = 0
            taxid_bp_dict = {}
            for line in report.readlines():
                column = line.split('\t')
                tax_id = column[2]
                read_len = int(column[3])
                total_bp += read_len
                if tax_id in taxid_bp_dict:
                    taxid_bp_dict[tax_id] += read_len
                else:
                    taxid_bp_dict[tax_id] = read_len

        with open(self._tool_inputs['TSV'][0].path) as tsv:

            expected = False
            allowed = False
            other = False
            for line in tsv.readlines():
                section = line.split('\t')
                taxa = section[5].strip()
                id = section[4]
                class_ = section[3].strip()
                percent = (taxid_bp_dict[id] / total_bp) * 100 if id in taxid_bp_dict else 0

                if class_ == self._parameters['level_of_depth'].value and (taxa == self._parameters['expected_species'].value or taxa in self._parameters['expected_species'].value):
                    expected = True
                    expected_spp = taxa
                    expected_percent = percent
                    continue

                if expected:
                    if class_ != 'S' and 'S' in class_:
                        expected_percent += percent
                        continue
                    else:
                        self._informs['expected'] = (expected_spp, round(expected_percent, 2))
                        expected = False

                if class_ == self._parameters['level_of_depth'].value and taxa in allowed_species:
                    allowed = True
                    allowed_spp = taxa
                    allowed_percent = percent
                    continue

                if allowed:
                    if class_ != 'S' and 'S' in class_:
                        allowed_percent += percent
                        continue
                    self._informs['allowed'].append((allowed_spp, round(allowed_percent, 2)))
                    allowed = False

                if not allowed and not expected and class_ == self._parameters['level_of_depth'].value:
                    other = True
                    other_spp = taxa
                    other_percent = percent
                    continue

                if other:
                    if class_ != 'S' and 'S' in class_:
                        other_percent += percent
                        continue
                    other = False
                    if other_percent > float(self._parameters['threshold_warn'].value):
                        if percent <= float(self._parameters['threshold_fail'].value):
                            self._informs['contaminants_warn'].append((other_spp, round(other_percent, 2),))
                        else:
                            self._informs['contaminants_fail'].append((other_spp, round(other_percent, 2),))
        if expected:
            self._informs['expected'] = (expected_spp, round(expected_percent, 2))
        elif allowed:
            self._informs['allowed'].append((allowed_spp, round(allowed_percent, 2)))
        elif other:
            if other_percent > float(self._parameters['threshold_warn'].value):
                if percent <= float(self._parameters['threshold_fail'].value):
                    self._informs['contaminants_warn'].append((other_spp, round(other_percent, 2),))
                else:
                    self._informs['contaminants_fail'].append((other_spp, round(other_percent, 2),))

        if 'expected' not in self._informs:
            logger.warning("No reads matching the expected species found! ")
            self._informs['expected'] = (self._parameters['expected_species'].value, 0)

        self._informs['contaminants_fail'].sort(key=lambda x: -x[1])
        self._informs['allowed'].sort(key=lambda x: -x[1])
        self._informs['contaminants_warn'].sort(key=lambda x: -x[1])
        self._informs['level_of_depth'] = self._parameters['level_of_depth'].value
        self._informs['threshold_fail'] = self._parameters['threshold_fail'].value
        self._informs['threshold_warn'] = self._parameters['threshold_warn'].value

    def _check_input(self) -> None:
        """
        Checks the tool input.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("TSV input is required.")
        super()._check_input()
