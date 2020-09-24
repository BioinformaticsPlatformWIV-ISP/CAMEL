import logging
from typing import List, Optional

import re

from camel.app.components.files.fastqutils import FastqUtils


class GalaxyUtils(object):
    """
    This class contains utility functions to work with Galaxy data sets.
    """

    @staticmethod
    def determine_sample_name_from_fq(fastq_names: List[str], is_pe: bool = True, default: Optional[str] = None) -> str:
        """
        Determines the sample name from the given command line arguments.
        :param fastq_names: PE FASTQ PE file names
        :param is_pe: If true, FASTQ files are paired-end
        :param default: Default value when the name cannot be parsed
        :return: Sample name
        """
        logging.debug(f"Determining sample name from: {', '.join(fastq_names)}")
        pattern = FastqUtils.PATTERN_FQ_PE if is_pe else FastqUtils.PATTERN_FQ_SE
        try:
            return FastqUtils.get_sample_name(fastq_names[0], pattern)
        except ValueError:
            logging.debug("Filename does not match any standard FASTQ format")

        # Trimmomatic output files
        m = re.search(r'.+ on {}'.format(pattern), fastq_names[0])
        if m:
            return m.group(1)

        # Raise error when default could not be set
        if default is None:
            raise ValueError(f"Sample name cannot be determined from: {', '.join(fastq_names)}")
        return default
