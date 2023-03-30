import logging
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class GmatsAlgorithm(Tool):

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('gMATS Algorithm', '0.1', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: concatenated tsv files
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(str(io.path) for io in self._tool_inputs['TSV']),
            '> /scratch/nagoeders/PyCharmProjects/RD_NG/camel_3.0/output_file.txt'
        ])
        self._execute_command()


if __name__ == '__main__':
    gmats = GmatsAlgorithm(Camel.get_instance())
    logging.info('test')
    gmats.add_input_files({'TSV': [
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-RRS16BD04259.tsv')),
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-S17BD02954.tsv')),
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-S17BD08805.tsv')),
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-S18BD04144.tsv')),
        ToolIOFile(Path('/testdata/camel/pipelines/gMATS/typing-bast-peptide-S18BD07986.tsv'))]
    })
    gmats.run()
