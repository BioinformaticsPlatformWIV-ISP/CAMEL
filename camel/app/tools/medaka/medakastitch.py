from camel.app.tools.medaka.medaka import Medaka


class MedakaStitch(Medaka):

    """
    Class for Medaka stitch function
    """

    def __init__(self, camel):
        """
        Initialize Medaka
        :param camel: Camel instance
        :return: None
        """
        super().__init__('medaka stitch', '1.7.3', camel)

        self._required_inputs = ['HDF', 'FASTA']
        self._output_type = 'FASTA'

    def _set_input(self):
        """
        Set the input specification
        :return: None
        """
        super(MedakaStitch, self)._set_input()

        hdf_file = self._tool_inputs['HDF'][0].path
        fasta_file = self._tool_inputs['FASTA'][0].path
        self._input_string = "{} {} ".format(hdf_file, fasta_file)
