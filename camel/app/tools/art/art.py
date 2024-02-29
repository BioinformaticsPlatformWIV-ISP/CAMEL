from pathlib import Path


from camel.app.tools.tool import Tool
from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile

class ART(Tool):

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool
        :param camel: CAMEL instance
        """
        super().__init__('ART', '2.5.8', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool
        :return: None
        """
        self.__build_command()
        self._execute_command()
        #TODO: self.__set_output()
        #TODO: self.__set_informs()

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

    def __build_command(self) -> None:
        """
        Builds the command. Only the art-illumina tool is implemented.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            '-ss', f"{self._parameters['sequencing_system'].value}",
            '-i', ' '.join(str(io.path) for io in self._tool_inputs['FASTA']),
            '-p',  # run a paired-end simulation
            '-l', f"{self._parameters['read_length'].value}",
            '-f', f"{self._parameters['fold_coverage'].value}",
            '-m', f"{self._parameters['mean_size'].value}",
            '-s', f"{self._parameters['standard_deviation'].value}",
            '-o', 'out'
            # TODO option '-na' is there to disable generating alignment files,
            # but this causes a warning, which camel somehow turns into an error.
        ])


if __name__ == '__main__':
    camel = Camel.get_instance()
    tool = ART(camel)
    tool.add_input_files({'FASTA': [ToolIOFile(Path('/db/refgenomes/Yersinia_enterocolitica/NC_GCA_02575835.1.fasta'))]})
    tool.run()