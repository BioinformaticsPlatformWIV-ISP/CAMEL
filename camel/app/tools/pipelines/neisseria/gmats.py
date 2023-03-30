import logging
from pathlib import Path

from camel.app.camel import Camel
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
        :return: None
        """
        print(f'{self.name} is running!')


if __name__ == '__main__':
    gmats = GmatsAlgorithm(Camel.get_instance())
    logging.info('test')
    gmats.add_input_files({'TSV' : [
        ToolIOFile(Path('typing-bast-peptide-RRS16BD04259.tsv')),
        ToolIOFile(Path('typing-bast-peptide-S17BD02954.tsv')),
        ToolIOFile(Path('typing-bast-peptide-S17BD08805.tsv')),
        ToolIOFile(Path('typing-bast-peptide-S18BD04144.tsv')),
        ToolIOFile(Path('typing-bast-peptide-S18BD07986.tsv'))]
    })
    gmats.run()