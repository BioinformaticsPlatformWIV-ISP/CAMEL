import abc
import logging
from pathlib import Path

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.tools.tool import Tool


class BcftoolsBase(Tool, metaclass=abc.ABCMeta):
    """
    Baseclass for bcftools.
    """

    def _get_output_key(self) -> str:
        """
        Returns the key for the output file.
        :return: Output file key
        """
        if 'output_type' not in self._parameters:
            return 'VCF'
        try:
            return {
                'b': 'BCF_GZ',
                'u': 'BCF',
                'z': 'VCF_GZ',
                'v': 'VCF'
            }[self._parameters['output_type'].value]
        except KeyError:
            raise InvalidParameterError(f"Invalid output format: {self._parameters['VCF'].value}")

    def _get_output_path(self) -> Path:
        """
        Returns the output filename.
        :return: Output filename
        """
        if 'output_filename' in self._parameters:
            return self.folder / FileSystemHelper.make_valid(self._parameters['output_filename'].value)
        logging.info("Output filename not set, reverting to default.")
        return self.folder / f"bcftools_out.{self._get_output_key().lower().replace('_', '.')}"

    def _check_command_output(self) -> None:
        """
        Checks if the command executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f'Error executing {self.name}: {self.stderr}')
