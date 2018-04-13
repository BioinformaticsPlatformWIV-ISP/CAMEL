import logging

import abc
import os

from app.components.vcf.vcfutils import VCFUtils
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class Filter(Tool, metaclass=abc.ABCMeta):
    """
    Base class for variant filters.
    """

    def __init__(self, *args):
        """
        Initializes this filter.
        :param args: Initialization arguments
        """
        super().__init__(*args)
        self._output_file = None

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if 'VCF_GZ' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No compressed VCF input found")
        super(Filter, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        nb_of_variants_pre = VCFUtils.count_variants(self._tool_inputs['VCF_GZ'][0].path)
        self._apply_filter()
        nb_of_variants_post = VCFUtils.count_variants(self.output_path)
        logging.info('{}/{} variants passed {} filtering'.format(nb_of_variants_post, nb_of_variants_pre, self._name))
        self._informs['variants_in'] = nb_of_variants_pre
        self._informs['variants_out'] = nb_of_variants_post
        self._tool_outputs['VCF_GZ'] = [ToolIOFile(self.output_path)]

    @abc.abstractmethod
    def _apply_filter(self):
        """
        Applies the filter to the input file.
        :return: None
        """
        pass

    @property
    def output_path(self):
        """
        Returns the path to the output file.
        :return: Path
        """
        return os.path.join(self._folder, self._parameters['output_filename'].value)
