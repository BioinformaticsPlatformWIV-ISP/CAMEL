from camel.app.tools.variantfiltering.basefilter import BaseFilter


class MappingQualityFilter(BaseFilter):
    """
    Filters variants based on mapping quality.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Variant Filter: Mapping Quality', '0.1', camel)

    def _apply_filter(self):
        """
        Applies the filtering on the variants.
        :return: None
        """
        self.__build_command()
        self._execute_command()

    @property
    def full_name(self) -> str:
        """
        Returns the full name for this filter.
        :return: Full name
        """
        return 'Mapping quality'

    def __build_command(self):
        """
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            "--exclude 'MQ<{}'".format(self._parameters['min_mapping_quality'].value),
            self._tool_inputs['VCF_GZ'][0].path,
            '--output-type z',
            '--output {}'.format(self.output_path)
        ])
