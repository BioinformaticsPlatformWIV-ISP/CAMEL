from camel.app.camel import Camel
from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4ApplyVQSR(GATK4):
    """
    =============================
    GATK ApplyVQSR 4.1.9.0
    =============================
    Apply a score cutoff to filter variants based on a recalibration table.

    Required inputs:
    ----------------
    'TXT_recal':            ToolIOFile object. The input recal file used by ApplyRecalibration
    'VCF':                  ToolIOFile object. One or more VCF files containing variants

    Output:
    -------
    'VCF'       The output filtered and recalibrated VCF file in which each variant is annotated with its VQSLOD value

    Mandatory parameters:
    ---------------------
    use_annotation:             The names of the annotations which should be used for calculations
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize GATKBaseRecalibrator tool.
        :param camel: Camel instance
        :return: None
        """
        super(GATK4BaseRecalibrator, self).__init__('gatk4 BaseRecalibrator', '4.1.9.0', camel)

        self._required_inputs = ['VCF', 'TXT_recal']
        self._output_type = 'TXT_RecalibrationTable'