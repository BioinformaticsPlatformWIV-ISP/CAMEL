import abc
from app.tools.tool import Tool


class Velvet(Tool, metaclass=abc.ABCMeta):

    """ Super class for Velvet related tools: velveth, VelvetOptimiser to handle inputs """

    FILEOPT_MAPPING = {
        'FASTA_REF': '-reference -fasta',
        'FASTA_long': '-long -fasta',
        'FASTQ_long': '-long -fastq',
        'FASTA_longPE': '-longPaired -fasta -separate',
        'FASTQ_longPE': '-longPaired -fastq -separate',
        'FASTA_SE': '-short -fasta',
        'FASTA_PE': '-shortPaired -fasta -separate',
        'FASTA_SE_2': '-short2 -fasta',
        'FASTA_SE_3': '-short3 -fasta',
        'FASTA_PE_2': '-shortPaired2 -fasta -separate',
        'FASTA_PE_3': '-shortPaired3 -fasta -separate',
        'FASTQ_SE': '-short -fastq',
        'FASTQ_PE': '-shortPaired -fastq -separate',
        'FASTQ_SE_2': '-short2 -fastq',
        'FASTQ_SE_3': '-short3 -fastq',
        'FASTQ_PE_2': '-shortPaired2 -fastq -separate',
        'FASTQ_PE_3': '-shortPaired3 -fastq -separate',
    }

    def __init__(self, tool, version, camel):
        """
        Initialize Velveth
        :param camel: Camel instance
        :param tool: the tool
        :param version: version of the tool
        :return: None
        """
        super(Velvet, self).__init__(tool, version, camel)
        self._input_string = None

    def _check_input(self):
        """
        Handle specific input for file handling
        :return: None
        """
        super(Velvet, self)._check_input()

        if not any(key in self._tool_inputs for key in ('FASTA_SE', 'FASTA_PE', 'FASTQ_SE', 'FASTQ_PE')):
            raise KeyError("No required PE or SE library input(s) found")

    def _set_input(self):
        """
        Set input
        :return: None
        """
        input_options = []

        for key, fileopt in Velvet.FILEOPT_MAPPING.items():
            if key in self._tool_inputs:
                if key.find('PE') > 0:
                    input_options.append("{} {}".format(fileopt, " ".join([f.path for f in self._tool_inputs[key]])))
                else:
                    input_options.append("{} {}".format(fileopt, self._tool_inputs[key][0].path))

        self._input_string = " ".join(input_options)
