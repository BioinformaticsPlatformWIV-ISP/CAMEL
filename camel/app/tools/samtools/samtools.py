from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.tools.toolpipeable import ToolPipeable


class Samtools(ToolPipeable):
    """
    Super class for samtools.
    """

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        super(Samtools, self)._execute_tool()

    def __init__(self, tool_name, version, camel):
        """
        Initialize a samtools tool.
        :param tool_name: Tool name
        :param version: Tool version
        :param camel: Camel instance
        :return: None
        """
        super().__init__(tool_name, version,  camel)

    def _check_stderr(self):
        """
        Validate if the program ran correctly by checking the standard error.
        :return: None
        """
        if any(keyword in self._command.stderr.lower() for keyword in ('aborted', 'error')):
            raise ToolExecutionError("{} failed: ''".format(self.name, self._command.stderr))
