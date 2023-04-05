from camel.app.tools.medaka.medaka import Medaka


class MedakaConsensus(Medaka):

    """
    Class for Medaka consensus function
    """

    def __init__(self, camel):
        """
        Initialize Medaka
        :param camel: Camel instance
        :return: None
        """
        super().__init__('medaka consensus', '1.7.3', camel)

        self._required_inputs = ['BAM']
        self._output_type = 'HDF'

    def _set_input(self):
        """
        Set the input specification
        :return: None
        """
        super(MedakaConsensus, self)._set_input()

        bam_file = self._tool_inputs['BAM'][0].path
        self._input_string = "{} ".format(bam_file)
