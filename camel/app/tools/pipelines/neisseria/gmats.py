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