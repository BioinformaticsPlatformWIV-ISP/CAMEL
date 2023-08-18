import re
from typing import Dict, Any, List

import bs4

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.tools.tool import Tool


class LREFinder(Tool):
    """
    Tool for detection of the 23S rRNA mutations, and the optrA, cfr, cfr(B) and poxtA genes, encoding linezolid
    resistance in Enterococci from whole genome sequences.
    Reference: https://cge.cbs.dtu.dk/services/LRE-Finder/
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('LRE-Finder', '20200812', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTQ_PE' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTQ_PE input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        output_path = self.folder / 'lrefinder_out.tsv'
        self._command.command = ' '.join([
            self._tool_command,
            '-ipe', *[str(f.path) for f in self._tool_inputs['FASTQ_PE']],
            '-t_db', '$LREFINDER_DB',
            '-o', str(output_path)
        ])
        self._execute_command()
        self._informs.update(self.__parse_html_output(self._command.stdout))

    # noinspection PyTypeChecker
    def __parse_html_output(self, html_code: str) -> Dict[str, Any]:
        """
        Parses the HTML output reported by the tool.
        :param html_code: HTML code
        :return: Parsed information
        """
        info = {}
        soup = bs4.BeautifulSoup(html_code, 'html.parser')
        info['species'] = soup.find('b', text='Species identified:').findNext('table').find('th').text
        info['genes'] = LREFinder.__parse_html_table(soup.find('b', text='Genes Identified:').findNext('table'))
        info['mutations'] = LREFinder.__parse_html_table(
            soup.find('b', text='Identified mutations in 23s:').findNext('table'))
        return info

    @staticmethod
    def __parse_html_table(html_table: bs4.ResultSet) -> List[Dict[str, Any]]:
        """
        Parses a HTML table.
        :param html_table: HTML table
        :return: Parsed information
        """
        # noinspection PyUnresolvedReferences
        rows = html_table.findAll('tr')
        header = [th.text.replace('_', ' ') for th in rows[0].findAll('th')]
        header = [re.sub(r'\[(.*)]', r'(\g<1>)', value) for value in header]
        table_data = [{
            col: td.text for col, td in zip(header, row.findAll('td'))} for row in rows[1:] if
            len(row.findAll('td')) > 0]

        # Convert numeric values to float
        for info_dict in table_data:
            for key in info_dict.keys():
                try:
                    info_dict[key] = float(info_dict[key])
                except ValueError:
                    continue
        return table_data

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if not self._command.returncode == 0:
            raise ToolExecutionError(self._command.stderr)
