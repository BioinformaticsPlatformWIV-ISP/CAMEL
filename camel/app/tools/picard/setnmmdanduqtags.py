from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.tools.picard.picard import Picard


class SetNmMdAndUqTags(Picard):

    """
    ==============================
    Picard SetNmMdAndUqTags 2.23.3
    ==============================
    Takes in a coordinate-sorted SAM or BAM and calculates the NM, MD, and UQ tags by comparing with the reference.

    Required inputs:
    ----------------
    ['BAM'|'SAM"]:      ToolIOFile object. A sorted SAM or BAM file.
    'FASTA_REF':        ToolIOFile object. FASTA file containing the reference genome.

    Output:
    -------
    'BAM'              ToolIOFile object. Fixed BAM file
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard SetNmMdAndUqTags ', '2.23.3', camel)

        self._function_name = 'SetNmMdAndUqTags '
        self._extra_inputs = ["FASTA_REF"]

    def _check_input(self) -> None:
        """
        Check input. Additionally, check whether required input is present
        :return: None
        """
        super(SetNmMdAndUqTags, self)._check_input()

        if 'FASTA_REF' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Picard SetNmMdAndUqTags: input file FASTA_REF is not defined")

    def _set_input(self) -> None:
        """
        Set input specification
        :return: None
        """
        super(SetNmMdAndUqTags, self)._set_input()

        self._input_string += f'R={self._tool_inputs["FASTA_REF"][0].path} '