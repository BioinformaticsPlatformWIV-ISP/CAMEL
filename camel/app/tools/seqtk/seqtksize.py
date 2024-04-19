from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.tools.tool import Tool


class SeqtkSize(Tool):

    """
    Reports the number sequences and bases in a FASTQ file.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize seqtk subsample
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Seqtk size', '1.4', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTQ' not in self._tool_inputs:
            raise InvalidInputSpecificationError(f'FASTQ input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Function to run seqtk subsample
        :return: None
        """
        self._informs['stats'] = []
        for path_fq in [x.path for x in self._tool_inputs['FASTQ']]:
            self._command.command = ' '.join([self._tool_command, str(path_fq)])
            self._execute_command()
            nb_reads, nb_bases = self.stdout.strip().split('\t')
            self.informs['stats'].append({
                'path': str(path_fq),
                'nb_of_bases': int(nb_bases),
                'nb_of_sequences': int(nb_reads)
            })
