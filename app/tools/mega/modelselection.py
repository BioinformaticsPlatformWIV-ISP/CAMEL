import logging
import os

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.invalidparametererror import InvalidParameterError
from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.mega.mltreeconstruction import MLTreeConstruction
from app.tools.tool import Tool


class ModelSelection(Tool):
    """
    Runs MEGA model selection.
    """

    DEFAULT_OUTPUT_NAME = 'model_selection'

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(ModelSelection, self).__init__('MEGA: Model Selection', '7.0.20', camel)

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No SNP Matrix FASTA input file found")
        super(ModelSelection, self)._check_input()

    def _check_parameters(self):
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
        if not os.path.isfile(self._parameters['config_file_template'].value):
            raise InvalidInputSpecificationError("Cannot read config file.")
        super(ModelSelection, self)._check_parameters()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        # self.__clear_output_files()
        self.__build_command()
        self._execute_command()
        self.__set_output()
        self.__analyze_output_file()

    def __clear_output_files(self):
        """
        Clears the output folder.
        :return: None
        """
        output_files = [os.path.join(self._folder, '{}.csv'.format(ModelSelection.DEFAULT_OUTPUT_NAME)),
                        os.path.join(self._folder, '{}_summary.txt'.format(ModelSelection.DEFAULT_OUTPUT_NAME))]
        for output_file in output_files:
            if os.path.isfile(output_file):
                os.remove(output_file)
                logging.debug("Removing '{}' from MEGA output folder".format(output_file))

    def __build_command(self):
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

    def __generate_config_file(self):
        """
        Generates the config file.
        :return: None
        """
        with open(self._parameters['config_file_template'].value) as handle:
            template = handle.read()

        config_file = os.path.join(self._folder, 'config.mao')
        with open(config_file, 'w') as handle:
            handle.write(template.format(
                branch_swap_filter=self._parameters['branch_swap_filter'].value,
                missing_data_treatment=self._parameters['missing_data_treatment'].value,
                site_coverage_cutoff=self._parameters['site_coverage_cutoff'].value if
                'site_coverage_cutoff' in self._parameters else 'Not Applicable'
            ))
        return config_file

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['CSV'] = [ToolIOFile(os.path.join(self._folder, '{}.csv'.format(
            ModelSelection.DEFAULT_OUTPUT_NAME)))]
        self._tool_outputs['TXT'] = [ToolIOFile(os.path.join(self._folder, '{}_summary.txt'.format(
            ModelSelection.DEFAULT_OUTPUT_NAME)))]

    def __analyze_output_file(self):
        """
        Analyzes the output file.
        :return: None
        """
        with open(self._tool_outputs['CSV'][0].path) as handle:
            self._informs['model'] = handle.readlines()[1].split(',')[0].split('+')[0]
            self._informs['model_full'] = MLTreeConstruction.SUBSTITUTION_MODELS[self._informs['model']]
            logging.info("Selected model: {}".format(self._informs['model']))

    def _check_command_output(self):
        """
        Checks the command output to see if the program executed correctly.
        :return: None
        """
        if 'MEGA-CC has logged the following error:' in self.stdout:
            raise ToolExecutionError("MEGA-CC failed to execute: {}".format(self.stdout.strip()))
