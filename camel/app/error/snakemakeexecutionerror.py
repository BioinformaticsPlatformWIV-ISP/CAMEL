from typing import Optional


class SnakemakeExecutionError(RuntimeError):
    """
    This error is raised when a Snakemake execution error occurs.
    """

    def __init__(self, stdout: str, stderr: str, failed_rule: Optional[str] = None) -> None:
        """
        This class is raised when a snakemake error occurs.
        :param stdout: Standard output
        :param stderr: Error output
        """
        super().__init__(f"Failed at rule: {failed_rule if failed_rule else 'n/a'}")
        self.stdout = stdout
        self.stderr = stderr
        self.failed_rule = failed_rule
