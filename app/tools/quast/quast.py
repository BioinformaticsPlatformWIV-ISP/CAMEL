import os

from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class Quast(Tool):
    """
    QUAST evaluates genome assemblies. QUAST works both with and without a reference genome. The tool accepts multiple
    assemblies, thus is suitable for comparison.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(Quast, self).__init__('quast', '4.4', camel)

    def _execute_tool(self):
        """
        Runs Quast
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA is required
        - FASTA_Ref, TSV_Gene, and TSV_Operon are optional
        - Only one input file allowed for FASTA_Ref, TSV_Gene, and TSV_Operon, multiple files allowed for FASTA
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise ValueError('No valid input key given for QUAST: {!r}'.format(self._tool_inputs))
        for key, values in self._tool_inputs.iteritems():
            if key not in ['FASTA', 'FASTA_Ref', 'TSV_Gene', 'TSV_Operon']:
                raise ValueError('Illegal input key given for QUAST: {!r}'.format(self._tool_inputs))
            if key in ['FASTA_Ref', 'TSV_Gene', 'TSV_Operon'] and len(values) > 1:
                raise ValueError('Too many input files given for QUAST: {!r}'.format(self._tool_inputs))
            for infile in values:
                if not infile.is_valid():
                    raise ValueError('Illegal input file given for QUAST: {!r}'.format(self._tool_inputs))

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command.
        :return: None
        """
        options_string = ' '.join(self._build_options() + ['-o {}'.format(self._folder)])
        input_string = self.__build_input_string()
        self._command.command = ' '.join([self._tool_command, options_string, input_string])

    def __build_input_string(self):
        """
        Creates the string with the input files
        :return: String with the input parameters
        """
        inputs = []
        if 'FASTA_Ref' in self._tool_inputs:
            inputs.append('-R {}'.format(self._tool_inputs['FASTA_Ref'][0].path))
        if 'TSV_Gene' in self._tool_inputs:
            inputs.append('-G {}'.format(self._tool_inputs['TSV_Gene'][0].path))
        if 'TSV_Operon' in self._tool_inputs:
            inputs.append('-O {}'.format(self._tool_inputs['TSV_Operon'][0].path))
        for item in self._tool_inputs['FASTA']:
            inputs.append(item.path)
        return ' '.join(inputs)

    def _check_command_output(self):
        """
        Checks if the command was executed successfully.
        :return: None
        """
        for line in self.stderr.splitlines():
            if 'ERROR' in line:
                if 'ERRORs: 0' not in line:
                    raise ToolExecutionError("Command execution failed (stderr: {}).".format(self.stderr))
        if self._command.returncode != 0:
            raise ToolExecutionError("Command execution failed (Exit code: {})".format(self._command.returncode))

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        output_keys = ['HTML', 'TEX', 'TSV', 'TXT']
        for key in output_keys:
            self._tool_outputs[key] = [ToolIOFile('{}.{}'.format(os.path.join(self._folder, 'report'), key.lower()))]
