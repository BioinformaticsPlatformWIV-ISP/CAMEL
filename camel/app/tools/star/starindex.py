from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.loggers import logger
from camel.app.tools.star.star import Star


class StarIndex(Star):
    """
    Index reference genome using STAR.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes STAR 2.7.11b indexing
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('STAR', '2.7.11b', camel)
        self._required_inputs = ['FASTA']
        self._input_string = "--runMode genomeGenerate"
        self._index_dir = ""

    def _set_input(self) -> None:
        """
        Set the input specification and the input string
        :return: None
        """
        self._input_string += " --genomeFastaFiles"
        for fasta in self._tool_inputs['FASTA']:
            self._input_string += f" {Path(str(fasta))}"

        if 'INDEX_DIR' not in self._tool_inputs:
            logger.warning("INDEX_DIR not specified; creating 'STAR_index' directory in same directory as FASTA input")
            self._index_dir = Path(str(self._tool_inputs['FASTA'][0])).parent / "STAR_index"
        else:
            self._index_dir = Path(str(self._tool_inputs['INDEX_DIR'][0]))
        self._index_dir.mkdir(parents=True, exist_ok=True)
        self._input_string += f" --genomeDir {self._index_dir}"

        if 'GTF' in self._tool_inputs:
            self._input_string += f" --sjdbGGTFfile {Path(str(self._tool_inputs['GTF'][0]))}"

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        self._tool_outputs['INDEX_DIR'] = [ToolIODirectory(Path(self._index_dir))]

    def _check_output(self) -> None:
        """
        Checks if the output is valid.
        :return: None
        """
        if not any(Path(str(self._tool_outputs['INDEX_DIR'][0])).iterdir()):
            raise IOError("'INDEX_DIR' is empty - index has not been created.")
        super()._check_output()
