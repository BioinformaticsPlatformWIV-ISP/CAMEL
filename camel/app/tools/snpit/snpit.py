from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class Snpit(Tool):
    """
    Whole genome SNP based identification of members of the Mycobacterium tuberculosis complex.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Snpit', '1.0', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'VCF' not in self._tool_inputs:
            raise InvalidInputSpecificationError('VCF input is required.')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_informs()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self.__generate_report_section())]

    def __build_command(self) -> None:
        """
        Builds the command line call.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            '--input {}'.format(self._tool_inputs['VCF'][0].path)
        ])

    def __set_informs(self) -> None:
        """
        Sets the informs for this tool.
        SNPit uses the following format string to report the output:
        %s\n%16s %16s %16s %.1f %%
        :return: None
        """
        lines = self._command.stdout.splitlines()
        if 'none identified' in lines[1]:
            for key in ('species', 'lineage', 'sublineage'):
                self._informs[key] = 'NA'
            self._informs['percent_matched'] = 0.0
            return
        try:
            self._informs['species'] = lines[1][0:16].strip()
            self._informs['lineage'] = lines[1][17:33].strip()
            self._informs['sublineage'] = lines[1][34:50].strip()
            self._informs['percent_matched'] = float(lines[1][51:-1])
        except IndexError:
            raise ToolExecutionError(f"Error executing {self.name}, cannot parse output: {lines}")

    def __generate_report_section(self) -> HtmlReportSection:
        """
        Creates a report section with the results of the tool.
        :return: Report section
        """
        section = HtmlReportSection('Snpit')
        table_data = [
            ['Species:', '<i>{}</i>'.format(self._informs['species']) if self._informs['species'] != '' else 'NA'],
            ['Lineage:', self._informs['lineage'] if self._informs['species'] != '' else 'NA'],
            ['Sublineage:', self._informs['sublineage'] if self._informs['sublineage'] != '' else 'NA'],
            ['Percent match:', '{:.2f}%'.format(self._informs['percent_matched'])]
        ]
        section.add_table(table_data, table_attributes=[('class', 'information')])
        return section

    def _check_command_output(self) -> None:
        """
        Checks if the command executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Error executing SNPit: {self._command.stderr}")
