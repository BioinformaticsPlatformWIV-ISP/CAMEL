import os

from app.io.tooliofile import ToolIOFile
from app.tools.bedtools.bedtools import Bedtools


class BedtoolsGenomecov(Bedtools):
    """
    Bedtools Genomecov func class
    """
    OUTPUT_FILE_BASENAME = 'bedtools_genomecov'

    def __init__(self, camel):
        """
        Initialize a samtools tool.
        :param camel: Camel instance
        :return: None
        """
        super(BedtoolsGenomecov, self).__init__('bedtools genomecov', '2.25.0', camel)
        self._output_filename = None
        self._required_inputs = ['BAM']

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__set_output()
        self.__build_command()
        self._execute_command()

    def __build_command(self):
        """
        Builds the command.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(self._build_options()),
            ' -ibam {}'.format(self._tool_inputs['BAM'][0].path),
            ' > {}'.format(self._output_filename)
        ])

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        if any(param in self._parameters.keys() for param in ['BedGraphWithZeroCoverage', 'BedGraph']):
            self._output_filename = BedtoolsGenomecov.OUTPUT_FILE_BASENAME + ".bed"
            self._tool_outputs['TXT_BED'] = [ToolIOFile(os.path.join(self._folder, self._output_filename))]

        if any(param in self._parameters.keys() for param in ['DepthWithZeroCoord', 'Depth']):
            self._output_filename = BedtoolsGenomecov.OUTPUT_FILE_BASENAME + ".tsv"
            self._tool_outputs['TSV'] = [ToolIOFile(os.path.join(self._folder, self._output_filename))]

    def _check_parameters(self):
        """
        Check the parameters
        :return: None
        """
        super(BedtoolsGenomecov, self)._check_parameters()

        # OUTPUT_FORMAT_OPTIONS are mutually exclusive, only ONE can be specified.
        OUTPUT_FORMAT_OPTIONS = ['BedGraphWithZeroCoverage', 'BedGraph', 'DepthWithZeroCoverage', 'Depth']

        params_excluded = []
        params_excluded += OUTPUT_FORMAT_OPTIONS
        output_opt_found = False

        for param in self._parameters.keys():
            if param in params_excluded:
                if not output_opt_found:
                    params_excluded.remove(param)
                    output_opt_found = True
                else:
                    raise ValueError("Only one output option should be specified.")

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        self._check_required_inputs()

        if len(self._tool_inputs['BAM']) != 1:
            raise ValueError("Exactly one BAM input file expected.")

        super(BedtoolsGenomecov, self)._check_input()
