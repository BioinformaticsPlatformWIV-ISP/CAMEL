import os

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.mega import TEMPLATE_MODEL_SELECT
from camel.app.tools.mega.mltreeconstruction import MLTreeConstruction
from camel.app.tools.tool import Tool


class ModelSelection(Tool):
    """
    Runs MEGA model selection.
    """

    DEFAULT_OUTPUT_NAME = 'model_selection'

    def __init__(self, camel: Camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('MEGA: Model Selection', '10.0.4', camel)

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No SNP Matrix FASTA input file found")
        super(ModelSelection, self)._check_input()

    def _check_parameters(self) -> None:
        """
        Checks if the parameters are valid.
        :return: None
        """
        if self._parameters['branch_swap_filter'].value not in (
                'None', 'Very Weak', 'Weak', 'Moderate', 'Strong', 'Very Strong'):
            raise InvalidParameterError("Branch swap filter parameter value is not valid.")
        if self._parameters['missing_data_treatment'].value not in (
                'Complete deletion', 'Use all sites', 'Partial deletion'):
            raise InvalidParameterError("Missing data treatment parameter value is not valid.")
        if self._parameters['missing_data_treatment'].value == 'Partial deletion':
            if 'site_coverage_cutoff' not in self._parameters:
                raise InvalidParameterError("No site coverage cutoff given for partial deletion")
            try:
                int(self._parameters['site_coverage_cutoff'].value)
            except ValueError:
                raise InvalidParameterError("Site coverage cutoff must be an integer")
        else:
            if 'site_coverage_cutoff' in self._parameters:
                raise InvalidParameterError("Site coverage cutoff is only applicable for 'Partial deletion'")
        if not os.path.isfile(TEMPLATE_MODEL_SELECT):
            raise InvalidInputSpecificationError("Cannot read config file.")
        super(ModelSelection, self)._check_parameters()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()
        self.__analyze_output_file()

    def __build_command(self) -> None:
        """
        Builds the command line call.
        :return: None
        """
        config_file = self.__generate_config_file()
        self._command.command = ' '.join([
            self._tool_command,
            '-d {}'.format(self._tool_inputs['FASTA'][0].path),
            '-a {}'.format(config_file),
            '-o {}'.format(ModelSelection.DEFAULT_OUTPUT_NAME)
        ])

    def __generate_config_file(self) -> str:
        """
        Generates the config file.
        :return: Path to config file
        """
        with open(TEMPLATE_MODEL_SELECT) as handle:
            template = handle.read()

        config_file = os.path.join(self._folder, 'config.mao')
        with open(config_file, 'w') as handle:
            handle.write(template.format(
                branch_swap_filter=self._parameters['branch_swap_filter'].value,
                missing_data_treatment=self._parameters['missing_data_treatment'].value,
                site_coverage_cutoff=self._parameters['site_coverage_cutoff'].value if
                'site_coverage_cutoff' in self._parameters else 'Not Applicable',
                threads=self._parameters['threads'].value
            ))
        return config_file

    def __set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['CSV'] = [ToolIOFile(os.path.join(self._folder, '{}.csv'.format(
            ModelSelection.DEFAULT_OUTPUT_NAME)))]
        self._tool_outputs['TXT'] = [ToolIOFile(os.path.join(self._folder, '{}_summary.txt'.format(
            ModelSelection.DEFAULT_OUTPUT_NAME)))]

    def __analyze_output_file(self) -> None:
        """
        Analyzes the output file.
        :return: None
        """
        with open(self._tool_outputs['CSV'][0].path) as handle:
            first_line = handle.readlines()[1]
            # The line has the following structure: K2+G+I, ...
            # The model is the first part (K2), +G means Gamma categories per site, +I means invariant sites
            complete_model = first_line.split(',')[0]
            self._informs['model'] = first_line.split(',')[0].split('+')[0]
            self._informs['model_full'] = MLTreeConstruction.SUBSTITUTION_MODELS[self._informs['model']]

            complete_rates = '+'.join(complete_model.split('+')[1:])
            self._informs['rates_among_sites'] = 'U' if complete_rates == '' else complete_rates
            self._informs['rates_among_sites_full'] = MLTreeConstruction.RATES_AMONG_SITES[self._informs[
                'rates_among_sites']]

    def _check_command_output(self) -> None:
        """
        Checks the command output to see if the program executed correctly.
        :return: None
        """
        if 'MEGA-CC has logged the following error:' in self.stdout:
            raise ToolExecutionError("MEGA-CC failed to execute: {}".format(self.stdout.strip()))
