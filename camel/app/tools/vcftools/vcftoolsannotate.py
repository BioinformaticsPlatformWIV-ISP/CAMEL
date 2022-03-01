from camel.app.camel import Camel
from camel.app.tools.vcftools.vcftools import VCFtools
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError


class VCFtoolsAnnotate(VCFtools):
    """
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize vcf-annotate tool.
        :param camel: Camel instance
        :return: None
        """
        super(VCFtools, self).__init__('VCFtools vcf-annotate', '0.1.16', camel)
        self._specific_parameters = ['output']

    def _build_command(self) -> None:
        """
        Build command of vcf-annotate function
        :return: None
        """
        build_options = self._build_options(excluded_parameters=self._specific_parameters)

        if 'VCF' in self._tool_inputs:
            input_string = f'cat {self._tool_inputs["VCF"][0].path}'
        elif 'VCF_GZ' in self._tool_inputs:
            input_string = f'gunzip -c {self._tool_inputs["VCF_GZ"][0].path}'
        else:
            raise InvalidInputSpecificationError("VCFtools vcf-annotate requires a VCF or VCF_GZ input file.")

        self._command.command = f'{input_string} | {self._tool_command} {" ".join(build_options)} | bgzip -c > {self._parameters["output"].value}'
