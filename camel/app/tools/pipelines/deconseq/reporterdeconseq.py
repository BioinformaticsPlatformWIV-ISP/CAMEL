from typing import List

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class ReporterDeconseq(Tool):
    """
    Tool to create HTML reports for deconseq.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize this tool.
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('Deconseq: reporter', '0.1', camel)
        # self.__sub_folder = 'read_trimming'
        self._report_section = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._report_section = HtmlReportSection('Deconseq', subtitle=self._input_informs['deconseq']['PE_FWD']['_name'])
        self.__add_deconseq_section()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section, False)]

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'deconseq' not in self._input_informs:
            raise InvalidInputSpecificationError("No trimming info found")
        super()._check_input()

    def __add_deconseq_section(self) -> None:
        """
        Adds the deconseq section.
        :return: None
        """
        self._report_section.add_header('Deconseq', 3)
        table_data = [self.__get_initial_read_counts()]
        for i, db in enumerate(self._input_informs['deconseq']['PE_FWD']['processed_dbs']):
            table_data.append(self.__get_db_removed_read_counts(db, i))
        table_data.append(self.__get_cleaned_read_counts())
        table_data.append(self.__get_processed_read_counts())
        self._report_section.add_table(table_data, column_names=['Metric', 'PE (fwd)', 'PE (rev)', 'SE'], table_attributes=[('class', 'data')])

    def __get_initial_read_counts(self) -> List[str]:
        """
        Returns the initial read counts for the PE and SE input files
        :return: List with initial read count values
        """
        read_counts = ['Initial reads']
        se_counts = 0
        for read_type in ['PE_FWD', 'PE_REV', 'SE_FWD', 'SE_REV']:
            if self._input_informs['deconseq'][read_type] is not None:
                if read_type.startswith('PE_'):
                    read_counts.append(ReporterDeconseq.__reformat_inform(str(self._input_informs['deconseq'][read_type]['initial_reads_count'])))
                else:
                    se_counts += self._input_informs['deconseq'][read_type]['initial_reads_count']
            elif read_type.startswith('PE_'):
                read_counts.append('-')
        read_counts.append(ReporterDeconseq.__reformat_inform(str(se_counts)))
        return read_counts

    def __get_db_removed_read_counts(self, db: str, index: int) -> List[str]:
        """
        Returns the read counts that were removed after running each db on the PE and SE input files
        :return: List with removed read count values
        """
        read_counts = [f'Removed {db} reads']
        se_counts = 0
        for read_type in ['PE_FWD', 'PE_REV', 'SE_FWD', 'SE_REV']:
            if self._input_informs['deconseq'][read_type] is not None:
                if read_type.startswith('PE_'):
                    read_counts.append(ReporterDeconseq.__reformat_inform(str(self._input_informs['deconseq'][read_type]['removed_reads_counts'][index])))
                else:
                    se_counts += self._input_informs['deconseq'][read_type]['removed_reads_counts'][index]
            elif read_type.startswith('PE_'):
                read_counts.append('-')
        read_counts.append(str(se_counts))
        return read_counts

    def __get_cleaned_read_counts(self):
        """
        Returns the read counts of the PE and SE input files that remain after processing all databases
        :return: List with cleaned read count values
        """
        read_counts = ['Cleaned reads (raw)']
        se_counts = 0
        for read_type in ['PE_FWD', 'PE_REV', 'SE_FWD', 'SE_REV']:
            if self._input_informs['deconseq'][read_type] is not None:
                if read_type.startswith('PE_'):
                    read_counts.append(ReporterDeconseq.__reformat_inform(str(self._input_informs['deconseq'][read_type]['final_reads_count'])))
                else:
                    se_counts += self._input_informs['deconseq'][read_type]['final_reads_count']
            elif read_type.startswith('PE_'):
                read_counts.append('-')
        read_counts.append(str(se_counts))
        return read_counts

    def __get_processed_read_counts(self):
        """
        Returns the read counts of the PE and SE input files that remain after processing all databases
        and after properly pairing all outputs again.
        :return: List with cleaned and completely processed read count values
        """
        read_counts = ['Finals reads (proper pairs)',
                       ReporterDeconseq.__reformat_inform(str(self._input_informs['deconseq']['combined']['remaining_pe_reads'])),
                       ReporterDeconseq.__reformat_inform(str(self._input_informs['deconseq']['combined']['remaining_pe_reads'])),
                       ReporterDeconseq.__reformat_inform(str(self._input_informs['deconseq']['combined']['remaining_se_reads']))]
        return read_counts

    @staticmethod
    def __reformat_inform(input_str: str) -> str:
        """
        This function is used to reformat an inform value to a more readable format.
        This function also works when the percentage is omitted.
        E.g. 5241241 (10.02%) -> 5.241.241 (10.02%)
        :param input_str: Input string
        :return: Reformatted inform
        """
        parts = input_str.split(' ')
        if len(parts) == 1:
            return '{:,}'.format(int(parts[0]))
        elif len(parts) == 2:
            return '{:,} {}'.format(int(parts[0]), parts[1])
        raise ValueError("Cannot parse: {}".format(input_str))
