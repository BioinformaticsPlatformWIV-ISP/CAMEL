from camel.app.tools.variantfiltering.filter import Filter


class SnpQualityFilter(Filter):
    """
    Filters variants based on SNP quality.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(SnpQualityFilter, self).__init__('Variant Filter: SNP Quality', '0.1', camel)

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
        return 'SNP quality'

    def __build_command(self):
        """
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            "--exclude 'QUAL<{}'".format(self._parameters['min_snp_quality'].value),
            self._tool_inputs['VCF_GZ'][0].path,
            '--output-type z',
            '--output {}'.format(self.output_path)
        ])
