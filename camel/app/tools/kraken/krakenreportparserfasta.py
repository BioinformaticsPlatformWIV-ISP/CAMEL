import pandas as pd

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class KrakenReportParserFasta(Tool):
    """
    Parses Kraken output reports for fasta input.
    """

    def __init__(self, camel) -> None:
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

        tsv_out = pd.read_table(self._tool_inputs['TSV_out'][0].path, header=None)
        total_bp = 0
        tax_id_bp_dict = {}
        for line in tsv_out.itertuples():
            tax_id = str(line[3])
            read_len = line[4]
            total_bp += read_len
            if tax_id in tax_id_bp_dict:
                tax_id_bp_dict[tax_id] += read_len
            else:
                tax_id_bp_dict[tax_id] = read_len

        tsv_report = pd.read_table(self._tool_inputs['TSV'][0].path, header=None)
        tsv_report.loc[len(tsv_report)] = ''
        expected = allowed = other = False
        expected_percent = allowed_percent = other_percent = 0
        expected_spp = allowed_spp = other_spp = None

        for section in tsv_report.itertuples():
            rank = section[4]
            tax_id = str(section[5])
            tax_name = section[6].strip()
            percent = (tax_id_bp_dict[tax_id] / total_bp) * 100 if tax_id in tax_id_bp_dict else 0
            percent = float(f'{percent:.2f}') if tax_id in tax_id_bp_dict else 0
            recheck_conditions = True

            while recheck_conditions:
                recheck_conditions = False
                if rank == self._parameters['level_of_depth'].value and (
                        tax_name in self._parameters['expected_species'].value):
                    expected = True
                    expected_spp = tax_name
                    expected_percent = percent

                elif expected:
                    if rank != self._parameters['level_of_depth'].value and (
                            self._parameters['level_of_depth'].value in rank or 'S' in rank):
                        expected_percent += percent
                    else:
                        self._informs['expected'] = (expected_spp, expected_percent)
                        expected = False
                        recheck_conditions = True

                elif rank == self._parameters['level_of_depth'].value and tax_name in allowed_species and not allowed:
                    allowed = True
                    allowed_spp = tax_name
                    allowed_percent = percent

                elif allowed:
                    if rank != self._parameters['level_of_depth'].value and (
                            self._parameters['level_of_depth'].value in rank or 'S' in rank):
                        allowed_percent += percent
                    else:
                        self._informs['allowed'].append((allowed_spp, allowed_percent))
                        allowed = False
                        recheck_conditions = True

                elif not expected and not allowed and not other and rank == self._parameters['level_of_depth'].value:
                    other = True
                    other_spp = tax_name
                    other_percent = percent

                elif other:
                    if rank != self._parameters['level_of_depth'].value and (
                            self._parameters['level_of_depth'].value in rank or 'S' in rank):
                        other_percent += percent
                        continue
                    else:
                        other = False
                        if other_percent > float(self._parameters['threshold_warn'].value):
                            if other_percent <= float(self._parameters['threshold_fail'].value):
                                self._informs['contaminants_warn'].append((other_spp, other_percent,))
                            else:
                                self._informs['contaminants_fail'].append((other_spp, other_percent,))
                        recheck_conditions = True

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
