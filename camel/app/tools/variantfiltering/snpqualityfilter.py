import logging
import os

from camel.app.components.vcf.vcfutils import VCFUtils
from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.variantfiltering.filter import Filter


class SnpQualityFilter(Filter):
    """
    Filters variants based on SNP quality.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(SnpQualityFilter, self).__init__('Variant Filter: SNP Quality', '0.1', camel)

    def _check_parameters(self):
        """
        Checks the command line parameters.
        :return: None
        """
        if 'min_snp_quality' not in self._parameters:
            raise InvalidParameterError("Parameter 'min_snp_quality' not found")
        super(SnpQualityFilter, self)._check_parameters()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        nb_of_variants_pre = VCFUtils.count_variants(self._tool_inputs['VCF_GZ'][0].path)
        self.__build_command()
        self._execute_command()
        output_file = os.path.join(self._folder, self._parameters['output_filename'].value)
        self._tool_outputs['VCF_GZ'] = [ToolIOFile(output_file)]
        nb_of_variants_post = VCFUtils.count_variants(output_file)
        logging.info('{}/{} variants passed snp quality filtering'.format(nb_of_variants_post, nb_of_variants_pre))
        self._informs['variants_in'] = nb_of_variants_pre
        self._informs['variants_out'] = nb_of_variants_post

    def __build_command(self):
        """
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            "--exclude 'QUAL<{}'".format(self._parameters['min_snp_quality'].value),
            self._tool_inputs['VCF_GZ'][0].path,
            '--output-type z',
            '--output {}'.format(self._parameters['output_filename'].value)
        ])
