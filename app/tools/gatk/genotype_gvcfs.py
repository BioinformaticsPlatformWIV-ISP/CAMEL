from app.tools.gatk.gatk import GATK


class GenotypeGVCFs(GATK):
    """
    Class for GATK GenotypeGVCFs function
    """

    def __init__(self, camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super(GenotypeGVCFs, self).__init__('gatk GenotypeGVCFs', '3.4.46', camel)

        self._function_name = 'GenotypeGVCFs'
        self._required_inputs = ['gVCF', 'FASTA_REF']
        self._output_type = 'VCF_MultipleSample'

    def _set_input(self):
        """
        Set the input specification
        :return: None
        """
        super(GenotypeGVCFs, self)._set_input()

        for f in self._tool_inputs['gVCF']:
            self._input_string += "--variant {} ".format(f.path)
