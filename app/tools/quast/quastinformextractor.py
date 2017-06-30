from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.tools.tool import Tool


class QuastInformExtractor(Tool):

    """
    Customized tool to parse QUAST TXT output to gather assembly QC information
    """
    COUNTING_STATIS = ['# N\'s per 100 kbp', '# mismatches per 100 kbp', '# indels per 100 kbp']
    CONTIG_STATIS = ['# contigs', 'N50', 'NG50', 'NGA50', 'Largest alignment', 'Largest contig']
    GENOME_STATIS = ['Total length', 'Total aligned length', 'Reference length',
                     'Unaligned length', 'Genome fraction (%)', 'Duplication ratio']
    GC_STATIS = ['GC (%)', 'Reference GC (%)']
    ABNORMAL_STATIS = ['# misassemblies', '# misassembled contigs',
                       'Misassembled contigs length', '# unaligned contigs']

    def __init__(self, camel):
        """
        Initialize QuastInformExtraction
        :param camel: Camel instance
        :return: None
        """
        super(QuastInformExtractor, self).__init__('Quast InformExtractor', '4.4', camel)

    def _execute_tool(self):
        """
        Entry point to run InterProScan
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
            raise InvalidInputSpecificationError('Quast TSV output required for QuastInformExtraction is missing.')

    def __set_informs(self):
        """
        Extract Quast QC information from input
        :return: None
        """
        qc_stats = {
            'counts': QuastInformExtractor.COUNTING_STATIS,
            'contig': QuastInformExtractor.CONTIG_STATIS,
            'genome': QuastInformExtractor.GENOME_STATIS,
            'gc': QuastInformExtractor.GC_STATIS,
            'abnormals': QuastInformExtractor.ABNORMAL_STATIS
        }
        for key in qc_stats.keys():
            self._informs[key] = {}

        with open(self._tool_inputs['TSV'][0].path, 'r') as inf:
            for line in [l.strip() for l in inf.readlines()]:
                if len(line.split("\t")) > 1:
                    stats_key, stats_value = line.split("\t")
                    for qc_sect, qc_sect_stats in qc_stats.items():
                        if stats_key in qc_sect_stats:
                            self._informs[qc_sect][stats_key] = stats_value
                            break
