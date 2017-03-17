import os

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class QuastInformExtracter(Tool):

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
        super(QuastInformExtracter, self).__init__('Quast InformExtracter', '4.4', camel)

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
        super(QuastInformExtracter, self)._check_input()
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError(
                'Quast TSV output required for QuastInformExtraction is missing. Tool_inputs : {!r}'.format(self._tool_inputs))

    def __set_informs(self):
        """
        Extract Quast QC information from input
        :return: None
        """
        self.informs['counts'] = {}
        self.informs['contig'] = {}
        self.informs['genome'] = {}
        self.informs['gc'] = {}
        self.informs['abnormals'] = {}
        with open(self._tool_inputs['TSV'][0].path, 'r') as inf:
            for l in inf.readlines():
                l = l.strip()
                if l.split("\t") > 1:
                    key, value = l.split("\t")
                    if key in QuastInformExtracter.COUNTING_STATIS:
                        self.informs['counts'][key] = value
                        continue
                    if key in QuastInformExtracter.CONTIG_STATIS:
                        self.informs['contig'][key] = value
                        continue
                    if key in QuastInformExtracter.GENOME_STATIS:
                        self.informs['genome'][key] = value
                        continue
                    if key in QuastInformExtracter.GC_STATIS:
                        self.informs['gc'][key] = value
                        continue
                    if key in QuastInformExtracter.ABNORMAL_STATIS:
                        self.informs['abnormals'][key] = value
                        continue
