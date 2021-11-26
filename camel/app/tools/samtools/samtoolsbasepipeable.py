import abc

from camel.app.camel import Camel
from camel.app.tools.samtools.samtoolsbase import SamtoolsBase
from camel.app.tools.toolpipeable import ToolPipeable


class SamtoolsBasePipeable(ToolPipeable, SamtoolsBase, metaclass=abc.ABCMeta):
    """
    Super class for pipeable samtools.
    """

    def __init__(self, tool_name: str, version: str, camel: Camel) -> None:
        """
        Initialize a samtools tool.
        :param tool_name: Tool name
        :param version: Tool version
        :param camel: Camel instance
        :return: None
        """
        super().__init__(tool_name, version,  camel)
