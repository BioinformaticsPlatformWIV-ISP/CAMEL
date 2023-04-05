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
            f"> {self._parameters['output_filename'].value}"
        ])
        self._execute_command()
        self._tool_outputs['TXT'] = [ToolIOFile(self.folder / self._parameters['output_filename'].value)]

    def _set_informs(self) -> None:
        """
        Collects the informs for this tool.
        """
        self._informs['nb_files'] = len(self._tool_inputs['TXT'])


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
    logging.info('test2')
    gmats.run()
    logging.info('test3')
    print(gmats.tool_outputs)
