import os

from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Blast(Tool):
    """
    Base class for BLAST tools.

    INPUT:
    - Query (FASTA): FASTA file
    - Subject (FASTA_Subject / DB_BLAST): Either a FASTA file or a BLAST database with the subject sequences
    """

    def __init__(self, tool_name, version, camel):
        """
        Initializes this tool.
        :param tool_name: Tool name
        :param version: Tool version
        :param camel: Camel instance
        """
        super(Blast, self).__init__(tool_name, version, camel)
        self.__subject_key = None

    def _check_input(self):
        """
        Checks whether the required input files are specified.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise ValueError('No FASTA input found')
        super(Blast, self)._check_input()

    def _execute_tool(self):
        """
        Runs Blast.
        :return: None
        """
        self.__subject_key = self.__get_subject_key()
        self.__output_key = self.__get_output_key()
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __get_subject_key(self):
        """
        Returns the key of the subject, this can be:
        - FASTA_Subject: FASTA file of the subject sequence
        - DB_BLAST: BLAST database created using makeblastdb.
        :return: Key
        """
        if all(key in self._tool_inputs for key in ['DB_BLAST', 'FASTA_Subject']):
            raise ValueError("Cannot use DB_BLAST and FASTA_Subject at the same time")
        elif 'DB_BLAST' in self._tool_inputs:
            return 'DB_BLAST'
        elif 'FASTA_Subject' in self._tool_inputs:
            return 'FASTA_Subject'
        else:
            raise ValueError("No subject (FASTA_Subject / DB_BLAST) found")

    def __get_output_key(self):
        """
        Returns the output key based on the output format.
        :return: Key
        """
        output_format = self._parameters['output_format'].value
        if output_format == '5':
            return 'XML'
        elif '6' in output_format or '7' in output_format:
            return 'TSV'
        elif output_format in ('8', '9', '11'):
            return 'ASN'
        elif '10' in output_format:
            return 'CSV'
        elif output_format == '12':
            return 'JSON'
        else:
            return 'TXT'

    def __build_command(self):
        """
        Builds the command line string.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            '-query {}'.format(self._tool_inputs['FASTA'][0].path),
            self.__get_subject_argument(),
            '-out {}'.format(self.__get_output_filename()),
            ' '.join(self._build_options(excluded_parameters=('output_filename',)))])

    def __get_subject_argument(self):
        """
        Returns the command line argument for the subject.
        :return: Command line argument
        """
        if self.__subject_key == 'FASTA_Subject':
            return '-subject {}'.format(self._tool_inputs['FASTA_Subject'][0].path)
        elif self.__subject_key == 'DB_BLAST':
            return '-db {}'.format(self._tool_inputs['DB_BLAST'][0].path)

    def __get_output_filename(self):
        """
        Returns the name of the output file.
        :return: Output name
        """
        if 'output_filename' in self._parameters:
            return self._parameters['output_filename'].value
        else:
            return self.__get_default_output_name()

    def __get_default_output_name(self):
        """
        Generates the default output name.
        :return: Output name
        """
        fasta_file_basename = os.path.splitext(os.path.basename(self._tool_inputs['FASTA'][0].path))[0]
        return '{}_{}.{}'.format(self._tool_command, fasta_file_basename, self.__get_output_key().lower())

    def __set_output(self):
        """
        Sets the output.
        :return: None
        """
        output_filename = os.path.join(self._folder, self.__get_output_filename())
        self._tool_outputs[self.__get_output_key()] = [ToolIOFile(output_filename)]

    def _check_command_output(self):
        """
        Checks the command output for errors.
        :return: None
        """
        if 'error' in self.stderr.lower() or self._command.returncode != 0:
            raise ToolExecutionError("Error executing {}: {}".format(self.name, self._command.stderr.strip()))
