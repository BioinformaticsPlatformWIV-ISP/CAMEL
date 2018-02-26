import logging
import os

from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class MLTreeConstruction(Tool):
    """
    Constructs maximum likelihood trees using MEGA.
    """

    DEFAULT_OUTPUT_NAME = 'tree_construction'

    SUBSTITUTION_MODELS = {
        'JC': 'Jukes-Cantor model',
        'K2': 'Kimura 2-parameter model',
        'T92': 'Tamura 3-parameter model',
        'HKY': 'Hasegawa-Kishino-Yano model',
        'TN93': 'Tamura-Nei model',
        'GTR': 'General Time Reversible model'
    }

    HEURISTIC_METHODS = {
        'NNI': 'Nearest-Neighbor-Interchange (NNI)',
        'SPR3': 'Subtree-Pruning-Regrafting - Fast (SPR level 3)',
        'SPR5': 'Subtree-Pruning-Regrafting - Extensive (SPR level 5)'
    }

    INITIAL_TREE = {
        'NJ_BioNJ': 'Make initial tree automatically (Default - NJ/BioNJ)',
        'Max_Parsimony': 'Make initial tree automatically (Maximum Parsimony)',
        'NJ': 'Make initial tree automatically (Neighbor Joining)',
        'BioNJ': 'Make initial tree automatically (BioNJ)'
    }

    RATES_AMONG_SITES = {
        'U': 'Uniform Rates',
        'G': 'Gamma Distributed (G)',
        'I': 'Has Invariant Sites (I)',
        'G+I': 'Gamma Distributed With Invariant Sites (G+I)'
    }

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(MLTreeConstruction, self).__init__('MEGA: ML Tree Construction', '7.0.20', camel)

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        super(MLTreeConstruction, self)._check_input()

    def _check_parameters(self):
        """
        Checks if the parameters are valid.
        :return: None
        """
        if self._parameters['test_of_phylogeny'].value not in ('None', 'Bootstrap method'):
            raise InvalidParameterError("Test of phylogeny must be either 'None' or 'Bootstrap method'")
        if self._parameters['test_of_phylogeny'].value == 'Bootstrap method':
            if 'bootstrap_replications' not in self._parameters:
                raise InvalidParameterError("Number of bootstrap replications has to be specified")
            try:
                int(self._parameters['bootstrap_replications'].value)
            except ValueError:
                raise InvalidParameterError("Number of bootstrap replications has to be an integer.")
        else:
            if 'bootstrap_replications' in self._parameters:
                raise InvalidParameterError("Number of bootstraps specified when no bootstraps are performed")
        if self._parameters['heuristic_method'].value not in MLTreeConstruction.HEURISTIC_METHODS:
            raise InvalidParameterError("Invalid heuristic method '{}' ({} supported)".format(
                self._parameters['heuristic_method'].value, ', '.join(MLTreeConstruction.HEURISTIC_METHODS.keys())))
        if self._parameters['initial_tree'].value not in MLTreeConstruction.INITIAL_TREE:
            raise InvalidParameterError("Invalid initial tree ({} supported)".format(
                ', '.join(MLTreeConstruction.INITIAL_TREE.keys())))
        if self._parameters['rates_among_sites'].value not in ('G', 'G+I') and 'gamma_categories' in self._parameters:
            raise InvalidParameterError("Gamma categories are only used when 'G' or 'G+I' rate models are used.")
        super(MLTreeConstruction, self)._check_parameters()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        # self.__clear_output_files()
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __clear_output_files(self):
        """
        Clears the output folder.
        :return: None
        """
        output_files = [os.path.join(self._folder, '{}.nwk'.format(MLTreeConstruction.DEFAULT_OUTPUT_NAME)),
                        os.path.join(self._folder, '{}_summary.txt'.format(MLTreeConstruction.DEFAULT_OUTPUT_NAME))]
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
            '-o {}'.format(MLTreeConstruction.DEFAULT_OUTPUT_NAME)
        ])

    def __get_parameter_value(self, name):
        """
        Returns the value of the parameter with the given name. Returns 'Not Applicable' if the parameter is not
        specified.
        :param name: Parameter name
        :return: Parameter value of 'Not Applicable'
        """
        if name in self._parameters:
            return self._parameters[name].value
        return 'Not Applicable'

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
                test_of_phylogeny=self.__get_parameter_value('test_of_phylogeny'),
                bootstrap_replications=self.__get_parameter_value('bootstrap_replications'),
                model=MLTreeConstruction.SUBSTITUTION_MODELS[self._parameters['model'].value],
                heuristic_method=MLTreeConstruction.HEURISTIC_METHODS[self._parameters['heuristic_method'].value],
                initial_tree=MLTreeConstruction.INITIAL_TREE[self._parameters['initial_tree'].value],
                rates_among_sites=MLTreeConstruction.RATES_AMONG_SITES[self._parameters['rates_among_sites'].value],
                branch_swap_filter=self.__get_parameter_value('branch_swap_filter'),
                gamma_categories=self.__get_parameter_value('gamma_categories'),
                missing_data_treatment=self.__get_parameter_value('missing_data_treatment'),
                site_coverage_cutoff=self.__get_parameter_value('site_coverage_cutoff'),
                threads=self.__get_parameter_value('threads')
            ))
        return config_file

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['NWK'] = [ToolIOFile(os.path.join(self._folder, '{}.nwk'.format(
            MLTreeConstruction.DEFAULT_OUTPUT_NAME)))]
        self._tool_outputs['TXT'] = [ToolIOFile(os.path.join(self._folder, '{}_summary.txt'.format(
            MLTreeConstruction.DEFAULT_OUTPUT_NAME)))]

    def _check_command_output(self):
        """
        Checks the command output to see if the tool executed correctly.
        :return: None
        """
        if 'error' in self.stdout.lower():
            raise ToolExecutionError("Problem generating tree using megacc: {}".format(
                '\n'.join(self.stdout.splitlines()[-5:])))
