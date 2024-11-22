from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.star.star import Star


class StarIndex(Star):
    """
    Index reference genome using STAR.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes STAR 2.7.11b indexing.
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('STAR', '2.7.11b', camel)
        self._required_inputs = ['FASTA']
        self._input_string = "--runMode genomeGenerate "
        self._index_dir = ""

    def _set_input(self) -> None:
        """
        Sets the input specification and the input string.
        :return: None
        """
        option_fasta = "--genomeFastaFiles"
        for fasta in self._tool_inputs['FASTA']:
            input_fasta = self._symlink_fasta(Path(str(fasta))) if 'symlink_input' in self._parameters else (Path(str(fasta)))
            option_fasta += f" {input_fasta}"

        option_gtf = ""
        if 'GTF' in self._tool_inputs:
            option_gtf = f"--sjdbGGTFfile {Path(str(self._tool_inputs['GTF'][0]))}"

        self._input_string += " ".join([option_fasta,
                                        option_gtf])

    def _set_output(self) -> None:
        """
        Sets the output specification.
        :return: None
        """
        index_dir = Path(str(self._tool_inputs['FASTA'][0])).parent / "GenomeDir"
        self._tool_outputs['INDEX_DIR'] = [ToolIODirectory(index_dir)]

    def _check_output(self) -> None:
        """
        Checks if the output is valid.
        :return: None
        """
        if not any(Path(str(self._tool_outputs['INDEX_DIR'][0])).iterdir()):
            raise IOError("INDEX_DIR is empty - index has not been created.")
        super()._check_output()

    def _symlink_fasta(self, fasta: Path) -> Path:
        """
        Creates a symlink for the fasta input. This avoids errors when there are no writing permissions on the directory of the input fasta.
        :param key: Input key
        :return: Path to symlink input
        """
        path_link = self._folder / fasta.name
        if not path_link.is_file():
            logger.info(f'Creating symlink for input file: {path_link}')
            path_link.symlink_to(fasta)
        return path_link


if __name__ == '__main__':
    star = StarIndex(Camel.get_instance())
    # star.add_input_files({'FASTA': [ToolIOFile(Path('/testdata/camel/star/S_20_721_prophage.fasta')),
    #                                 ToolIOFile(Path('/testdata/camel/star/pBAD33.fasta'))],
    #                       'INDEX_DIR': [ToolIODirectory(Path('/scratch/grdeclercq/star/index'))],
    #                       })
    star.add_input_files({'FASTA': [ToolIOFile(Path('/scratch/grdeclercq/star/test/S_20_721_prophage.fasta'))]})
    star.update_parameters(SA_index=7)
    star.run(Path('/scratch/grdeclercq/star/test'))
