import re
from pathlib import Path

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.sequencetyping.sequencetypingsrst2hit import SequenceTypingSRST2Hit
from camel.app.components.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class SRST2AlleleDetector(Tool):
    """
    This tool detects an allele using SRST2.
    """

    __COLUMN_INDICES = {'allele_id': 2, 'mismatches': 3, 'uncertainty': 4, 'depth': 5}

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance.
        """
        super().__init__("SRST2: Allele detector", "0.2.0", camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        if 'FASTQ_PE' in self._tool_inputs:
            input_str = '--input_pe {}'.format(' '.join([str(x.path) for x in self._tool_inputs['FASTQ_PE']]))
        else:
            input_str = '--input_se {}'.format(self._tool_inputs['FASTQ_SE'][0].path)
        output_prefix = self._folder / f"_detection_{FileSystemHelper.make_valid(self._input_informs['locus']['name'])}"
        self._command.command = '{} {} --log --mlst_db {} --mlst_delimiter "{}" --output {} {}'.format(
            self._tool_command,
            input_str,
            self._tool_inputs['FASTA'][0].path,
            SequenceTypingUtils.determine_delimiter(self._tool_inputs['FASTA'][0].path),
            output_prefix,
            ' '.join(self._build_options())
        )
        self._execute_command()

        # Collect the output
        try:
            output_file = next(self.folder.glob('*__results.txt'))
        except StopIteration:
            raise ToolExecutionError("No output file generate by SRST2")
        hit = self.__extract_hit(self._input_informs['locus']['name'], output_file)
        self._tool_outputs['VAL_Hit'] = [ToolIOValue(hit)]

    def _check_input(self) -> None:
        """
        Checks if the required input is specified.
        :return: None
        """
        if ('FASTQ_PE' not in self._tool_inputs) and ('FASTQ_SE' not in self._tool_inputs):
            raise InvalidInputSpecificationError("FASTQ input (PE / SE) is required")
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Allele FASTA file is required")
        if 'locus' not in self._input_informs:
            raise InvalidInputSpecificationError("Locus information is required")
        super()._check_input()

    def _check_command_output(self) -> None:
        """
        Checks the command output.
        :return: None
        """
        if 'Could not determine forward/reverse read status' in self._command.stderr:
            raise ToolExecutionError("Invalid names for the FASTQ files")
        elif self._command.returncode != 0:
            raise ToolExecutionError("Error executing SRST2")

    def __extract_hit(self, locus_name: str, srst2_output: Path) -> SequenceTypingSRST2Hit:
        """
        Extracts a SRST2TypingHit from the SRST2 output.
        :param locus_name: Locus name
        :param srst2_output: SRST2 output file
        :return: SRST2 typing hit
        """
        with srst2_output.open() as handle:
            lines = handle.readlines()
            if len(lines) != 2:
                return SequenceTypingSRST2Hit.create_empty_hit(locus_name)
            parts = lines[1].strip().split('\t')
            return SequenceTypingSRST2Hit(
                locus_name,
                re.sub('[*?]', '', parts[self.__COLUMN_INDICES['allele_id']]),
                self.__clean_srst2_column(parts[self.__COLUMN_INDICES['mismatches']]),
                self.__clean_srst2_column(parts[self.__COLUMN_INDICES['uncertainty']]),
                float(parts[self.__COLUMN_INDICES['depth']]) if parts[self.__COLUMN_INDICES['depth']] != '-' else None
            )

    @staticmethod
    def __clean_srst2_column(input_column: str) -> str:
        """
        Cleans the SRST2 output.
        The mismatches and uncertainty column typically contain the gene name which have to be removed for the hit.
        Example:
            folP_3\\edge2.0 -> edge2.0
        :param input_column: Input column
        :return: Cleaned column
        """
        return re.sub('^.+/', '', input_column)
