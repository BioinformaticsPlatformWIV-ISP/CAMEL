from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.tools.tool import Tool


class BcftoolsIndexStats(Tool):
    """
    Returns stats based on existing index files.
    Stats:
    - Number of variants per reference sequence ('variants_per_reference')
    - Total number of variants ('total_variants')
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance.
        """
        super(BcftoolsIndexStats, self).__init__('bcftools index stats', '1.3.1', camel)

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('BCF', 'VCF_GZ')):
            raise InvalidInputSpecificationError("No input file found")
        if len(self._tool_inputs) != 1:
            raise InvalidInputSpecificationError("Only one type of input is supported (VCF_GZ or BCF)")
        super(BcftoolsIndexStats, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self._input_key = self._tool_inputs.keys()[0]
        self._command.command = '{} {}'.format(
            self._tool_command,
            self._tool_inputs[self._input_key][0].path
        )
        self._execute_command()
        self.__analyze_stdout()

    def __analyze_stdout(self):
        """
        Analyzes the standard output to report the SNP statistics.
        :return: None
        """
        self._informs['variants_per_reference'] = {}
        for line in self.stdout.splitlines():
            reference, size, snps = line.split('\t')
            self._informs['variants_per_reference'][reference] = int(snps)
        self._informs['total_variants'] = sum(self._informs['variants_per_reference'][ref] for ref in
                                              self._informs['variants_per_reference'].keys())
