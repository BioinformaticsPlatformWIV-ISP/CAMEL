import math
from pathlib import Path
from typing import Tuple, List

import numpy as np
import pandas as pd

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class PrimerProbeCheckReporterCSV(Tool):
    """
    Tool to create CSV files with information on the number of matches between primer/probe sequences and a fasta file.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize this tool.
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('Primer Probe Check: reporter', '0.1', camel)
        self._param = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._check_input()
        self._retrieve_parameters()
        self.__set_output(self._parameters['output_stat'].value, self._parameters['output_full_report'].value)

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'CSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No CSV files found")
        super()._check_input()

    def _retrieve_parameters(self) -> None:
        """
        Retrieve required parameters
        :return: None
        """
        self.primer = str(self._parameters['primer'].value)
        self.fasta_primer_name = str(self._parameters['fasta_primer_name'].value)
        self.end_mismatch = int(self._parameters['end_mismatch'].value)
        self.perc_mismatch = float(self._parameters['perc_mismatch'].value)

    @staticmethod
    def split_and_generate_names(data_list, len_excel, prefix) -> Tuple[List[pd.DataFrame], List[str]]:
        """
        There is a maximum number of rows to an Excel sheet, therefore, the dataframe is split into multiple sheets.
        Moreover, a unique name is generated with a given prefix.
        :param data_list: Dataframe that has to be split
        :param len_excel: Number of rows in Excel sheet
        :param prefix:  Prefix for the name of the sheet
        :return: Tuple containing a list of data chunks and a list of generated names.
        """
        df_list = np.array_split(data_list, math.ceil(data_list.shape[0] / len_excel))
        names = [f'{prefix}_{i}' for i in range(1, len(df_list) + 1)]
        return df_list, names

    @staticmethod
    def write_data_to_excel(data_list, writer, sheet_names) -> None:
        """
        Writes data to an Excel file using pandas.
        :param data_list: Dataframe that has to be split
        :param writer: Output filename
        :param sheet_names: Names of sheet(s)
        :return: None
        """
        for i, data in enumerate(data_list):
            data.to_excel(writer, sheet_name=sheet_names[i], index=False)

    def __csv_to_xlsx_full_report(self, outfile_name_stat: str, outfile_name_full_report: str) -> None:
        """
        Downstream analysis of CSV (combined) files to get the summary stats and the full report.
        :param outfile_name_stat: String of path output for summary statistics
        :param outfile_name_full_report: String of path output for full report
        :return: None
        """
        #Read CSV file
        df_csv = pd.read_csv(str(self._tool_inputs['CSV'][0]))
        # Determine mismatch_end_final based on fasta_primer_name
        mismatch_end_final = 0 if 'PR' in self.fasta_primer_name else self.end_mismatch

        # Remove duplicate rows based on 'id' column
        df_filtered = df_csv.drop_duplicates(subset='id')
        # Remove NaN rows from df_csv and remove duplicate rows from non-NaN rows
        df_non_nan_rows = df_csv[~df_csv['matched'].isna()]
        # Remove duplicate rows based on 'id' column
        df_non_nan_rows = df_non_nan_rows.drop_duplicates(subset='id')

        # Create a DataFrame 'stat_df' with summary statistics
        stat_df = pd.DataFrame({'Sequence_primer': [self.fasta_primer_name],
                                'Sequence': [str(self.primer)],
                                'Max_Percentage': [self.perc_mismatch],
                                'Number_Mismatch_3_end': [mismatch_end_final],
                                'Inclusivity': [df_non_nan_rows['matched'].count() / df_filtered.shape[0]],
                                'FN': [df_filtered.shape[0] - df_non_nan_rows['matched'].count()],
                                'Total': [df_filtered.shape[0]]})

        with pd.ExcelWriter(str(outfile_name_stat)) as writer:
            stat_df.to_excel(writer, sheet_name='Summary Results', index=False)

        # Separate matched and non-matched records
        match_list = df_csv[df_csv['matched'].notnull()]
        no_match_list = df_csv[df_csv['matched'].isnull()]

        # Define chunk size because Excel has a limit for its rows
        len_excel = 750000

        with pd.ExcelWriter(str(outfile_name_full_report)) as writer:
            # Write summary results
            stat_df.to_excel(writer, sheet_name='Summary Results', index=False)

            # Write match_list chunks
            if len(match_list) > 0:
                df_list_match, names_match = self.split_and_generate_names(match_list, len_excel, 'Match_list')
                self.write_data_to_excel(df_list_match, writer, names_match)

            # Write no_match_list chunks
            if len(no_match_list) > 0:
                df_list_no_match, names_nomatch = self.split_and_generate_names(no_match_list, len_excel,
                                                                                'No_Match_list')
                self.write_data_to_excel(df_list_no_match, writer, names_nomatch)

    def __set_output(self, output_stat: str, output_full_report: str) -> None:
        """
        Sets the tool output.
        :return: None
        """
        self._tool_outputs['CSV_STAT'] = [ToolIOFile(Path(output_stat))]
        self._tool_outputs['CSV_FP'] = [ToolIOFile(Path(output_full_report))]
        self.__csv_to_xlsx_full_report(output_stat, output_full_report)
