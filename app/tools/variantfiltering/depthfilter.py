import logging

import os

from app.components.vcf import vcfutils
from app.error.invalidparametererror import InvalidParameterError
from app.io.tooliofile import ToolIOFile
from app.tools.variantfiltering.filter import Filter


class DepthFilter(Filter):
    """
    Filters variants based on absolute depth, forward depth and reverse depth.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(DepthFilter, self).__init__('Variant Filter: Depth', '0.1', camel)

    def _check_parameters(self):
        """
        Checks the command line parameters.
        :return: None
        """
        if not any([param in self._parameters for param in ('min_depth', 'min_relative_depth', 'min_reverse_depth')]):
            raise InvalidParameterError('No filtering parameter found')
        super(DepthFilter, self)._check_parameters()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        nb_of_variants_pre = vcfutils.count_variants(self._tool_inputs['VCF_GZ'][0].path)
        self.__build_command()
        self._execute_command()
        output_file = os.path.join(self._folder, self._parameters['output_filename'].value)
        self._tool_outputs['VCF_GZ'] = [ToolIOFile(output_file)]
        nb_of_variants_post = vcfutils.count_variants(output_file)
        logging.info('{}/{} variants passed depth filtering'.format(nb_of_variants_post, nb_of_variants_pre))
        self._informs['variants_in'] = nb_of_variants_pre
        self._informs['variants_out'] = nb_of_variants_post

    def __create_exclude_string(self):
        """
        Creates the exclude string.
        :return: Exclude string
        """
        parts = []
        if 'min_depth' in self._parameters:
            parts.append('DP<{}'.format(self._parameters['min_depth'].value))
        if 'min_forward_depth' in self._parameters:
            parts.append('DP4[0]+DP4[2]<{}'.format(self._parameters['min_forward_depth'].value))
        if 'min_reverse_depth' in self._parameters:
            parts.append('DP4[1]+DP4[3]<{}'.format(self._parameters['min_reverse_depth'].value))
        return ' || '.join(parts)

    def __build_command(self):
        """
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            "--exclude '{}'".format(self.__create_exclude_string()),
            self._tool_inputs['VCF_GZ'][0].path,
            '--output-type z',
            '--output {}'.format(self._parameters['output_filename'].value)
        ])
