import re
from pathlib import Path

from camel.app.tools.tool import Tool
from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.error.toolexecutionerror import ToolExecutionError



class ART(Tool):
    """
    A simulation tool to generate synthetic next-generation sequencing reads.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('ART', '2.5.8', camel)
        #TODO: support single-end read simulation mode as well?

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Build the ART command and execute it
        self.__build_command()
        self._execute_command()
        # Compress the two output files
        FileSystemHelper.gzip_file(self.__get_output_path('1'), Path(f"{self.__get_output_path('1')}.gz"))
        FileSystemHelper.gzip_file(self.__get_output_path('2'), Path(f"{self.__get_output_path('2')}.gz"))
        #TODO: delete .fq files?
        # Set camel output and informs
        self.__set_output()
        self.__set_informs()
        if not self._informs.get('succeed', 'False'):
            raise ToolExecutionError("Error running ART")

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if 'FASTA' in self._tool_inputs:
            if len(self._tool_inputs['FASTA']) != 1:
                raise ValueError("FASTA input requires exactly 1 file.")
        else:
            raise ValueError("No FASTA input found")
        super(ART, self)._check_input()

    def __build_command(self) -> None:
        """
        Builds the command.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            '-ss', f"{self._parameters['sequencing_system'].value}",
            '-i', ' '.join(str(io.path) for io in self._tool_inputs['FASTA']),
            '-p', # Run a paired-end simulation
            '-l', f"{self._parameters['read_length'].value}",
            '-f', f"{self._parameters['fold_coverage'].value}",
            '-m', f"{self._parameters['mean_size'].value}",
            '-s', f"{self._parameters['standard_deviation'].value}",
            '-o', f"{self._parameters['out'].value}",
            '-na' # Do not generate an ALN alignment file
        ])

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if re.match("Warning: your simulation will not output any ALN", self.stderr):
            # The warning can be ignored because we don't need ALN files
            pass
        else:
            super(ART, self)._check_command_output()

    def __get_output_path(self, suffix: str) -> Path:
        """
        Returns the path for the output file with the given suffix.
        """
        basename = self._parameters['out'].value
        return self.folder / f"{FileSystemHelper.make_valid(basename)}{suffix}.fq"

    def __set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['FASTQ'] = [
            ToolIOFile(Path(f"{self.__get_output_path('1')}.gz")), ToolIOFile(Path(f"{self.__get_output_path('2')}.gz"))]

    def __set_informs(self) -> None:
        """
        Adds to the informs whether execution succeeded.
        :return: None
        """
        #TODO: is there a better way to do this? --> not necessary
        for line in self.stdout.splitlines():
            if line.strip() == 'Output files':
                self._informs['succeed'] = True

if __name__ == '__main__':
    camel = Camel.get_instance()
    tool = ART(camel)
    tool.add_input_files({'FASTA': [ToolIOFile(Path('/db/refgenomes/Yersinia_enterocolitica/NC_GCA_02575835.1.fasta'))]})
    tool.run()
    print("OUTPUT: ",tool._tool_outputs)
    print("INFORMS: ", tool._informs)