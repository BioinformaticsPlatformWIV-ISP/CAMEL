from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.tools.tool import Tool


class QuastInformExtractor(Tool):

    """
    Customized tool to parse QUAST TXT output to gather assembly QC information
    """

    KEY_MAPPING = {
        'counts': (
            '# N\'s per 100 kbp',
            '# mismatches per 100 kbp',
            '# indels per 100 kbp'),
        'contig': (
            '# contigs',
            '# contigs (>= 0 bp)',
            '# contigs (>= 1000 bp)',
            'N50',
            'NG50',
            'NGA50',
            'Largest alignment',
            'Largest contig'),
        'genome': (
            'Total length',
            'Total length (>= 0bp)',
            'Total length (>= 1000bp)',
            'Total aligned length',
            'Reference length',
            'Unaligned length',
            'Genome fraction (%)',
            'Duplication ratio'),
        'gc': (
            'GC (%)',
            'Reference GC (%)'),
        'abnormal': (
            '# misassemblies',
            '# misassembled contigs',
            'Misassembled contigs length',
            '# unaligned contigs')
    }

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Quast InformExtractor', '4.4', camel)

    def _execute_tool(self):
        """
        Execute this tool.
        :return: None
        """
        self.__set_informs()

    def _check_input(self):
        """
        Checks whether required quast TSV input is available
        :return: None
        """
        super(QuastInformExtractor, self)._check_input()
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("TSV input file is required.")

    def __set_informs(self):
        """
        Extract Quast QC information from input file.
        :return: None
        """
        self._informs = {k: {} for k in QuastInformExtractor.KEY_MAPPING.keys()}
        with open(self._tool_inputs['TSV'][0].path, 'r') as handle:
            for key, value in [l.strip().split('\t') for l in handle.readlines()]:
                for mapping_name, mapping_keys in QuastInformExtractor.KEY_MAPPING.items():
                    if key in mapping_keys:
                        self._informs[mapping_name][key] = value
                        break
