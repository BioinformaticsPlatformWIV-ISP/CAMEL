import json
import re
from typing import Tuple, Dict


class GeneDetectionUtils(object):
    """
    This class contains utility functions for the gene detection workflow.
    """

    @staticmethod
    def parse_header(header: str) -> Tuple[str, Dict]:
        """
        Parses a gene detection header. The format is:
        >{sequence id} {metadata in JSON format}
        :param header: Complete header
        :return: sequence id, metadata
        """
        m = re.match('^(.*) ({.*})$', header)
        if not m:
            raise ValueError("Invalid header: {}".format(header))
        metadata = json.loads(m.group(2))
        return m.group(1), metadata
