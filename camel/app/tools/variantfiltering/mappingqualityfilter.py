from app.tools.variantfiltering.filter import Filter


class MappingQualityFilter(Filter):
    """
    Filters variants based on mapping quality.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(MappingQualityFilter, self).__init__('Variant Filter: Mapping Quality', '0.1', camel)

    def _apply_filter(self):
        """
        Applies the filtering on the variants.
        :return: None
        """
        self.__build_command()
        self._execute_command()

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
