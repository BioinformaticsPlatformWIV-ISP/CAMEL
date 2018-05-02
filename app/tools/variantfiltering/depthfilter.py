from app.error.invalidparametererror import InvalidParameterError
from app.tools.variantfiltering.filter import Filter


class DepthFilter(Filter):
    """
    Filters variants based on absolute depth, forward depth and reverse depth.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(DepthFilter, self).__init__('Variant Filter: Depth', '0.1', camel)

    def _check_parameters(self):
        """
        Checks the command line parameters.
        :return: None
        """
        if not any([param in self._parameters for param in ('min_depth', 'min_relative_depth', 'min_reverse_depth')]):
            raise InvalidParameterError('No filtering parameter found')
        super(DepthFilter, self)._check_parameters()

    def _apply_filter(self):
        """
        Applies the filtering on the variants.
        :return: None
        """
        self.__build_command()
        self._execute_command()

    def __create_exclude_string(self):
        """
        Creates the exclude string.
        :return: Exclude string
        """
        parts = []
        if 'min_depth' in self._parameters:
            parts.append('DP<{}'.format(self._parameters['min_depth'].value))
        if 'min_forward_depth' in self._parameters:
            parts.append('DP4[0]+DP4[2]<{}'.format(self._parameters['min_forward_depth'].value))
        if 'min_reverse_depth' in self._parameters:
            parts.append('DP4[1]+DP4[3]<{}'.format(self._parameters['min_reverse_depth'].value))
        return ' || '.join(parts)

    def __build_command(self):
        """
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            "--exclude '{}'".format(self.__create_exclude_string()),
            self._tool_inputs['VCF_GZ'][0].path,
            '--output-type z',
            '--output {}'.format(self.output_path)
        ])
