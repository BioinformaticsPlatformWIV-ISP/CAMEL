import os
from abc import ABCMeta

import logging

from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class BlastLoop(Tool, metaclass=ABCMeta):
    """
    Tool that loops over input files to run BLAST multiple times.
    """

    def __init__(self, tool_name, version, camel):
        """
        Initializes this tool.
        :param tool_name: Tool name
        :param version: Tool version
        :param camel: Camel instance
        """
        super(BlastLoop, self).__init__(tool_name, version, camel)
        self.__subject_key = None

    def _check_input(self):
        """
        Checks whether the required input files are specified.
        :return: None
        """
        if len(self._tool_inputs) <= 1:
            pass
        elif 'FASTA' not in self._tool_inputs:
            raise ValueError('No FASTA input found')
        super(BlastLoop, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        if len(self._tool_inputs) <= 1:
            logging.info("Not enough inputs to run {}".format(self.name))
            return
        self.__subject_key = self.__get_subject_key()
        self._tool_outputs[self.__get_output_key()] = []

        if len(self._tool_inputs['FASTA']) > 1 and len(self._tool_inputs[self.__subject_key]) > 1:
            raise ValueError("Cannot loop over both inputs")
        elif len(self._tool_inputs['FASTA']) > 1:
            subject = self._tool_inputs[self.__subject_key][0]
            for query in self._tool_inputs['FASTA']:
                self.__run_blast(query, subject)

        elif len(self._tool_inputs[self.__subject_key]) > 1:
            query = self._tool_inputs['FASTA'][0]
            for subject in self._tool_inputs[self.__subject_key]:
                self.__run_blast(query, subject)
        else:
            self.__run_blast(self._tool_inputs['FASTA'][0], self._tool_inputs[self.__subject_key][0])

    def __run_blast(self, query, subject):
        """
        Runs blast.
        :return: None
        """
        self._command.command = self.__build_command(query, subject)
        self._execute_command()
        self._tool_outputs[self.__get_output_key()].append(
            ToolIOFile(os.path.join(self._folder, self.__get_output_name(query, subject))))

    def __build_command(self, query, subject):
        """
        Builds the command line command.
        :param query: Query
        :param subject: Subject
        :return: Command line string
        """
        return ' '.join([self._tool_command,
                         '-query {}'.format(query.path),
                         self.__get_subject_argument(subject.path),
                         '-out {}'.format(self.__get_output_name(query, subject)),
                         ' '.join(self._build_options())])

    def __get_subject_argument(self, subject_file):
        """
        Returns the command line argument for the subject.
        :param subject_file: Subject file
        :return: Command line argument
        """
        if self.__subject_key == 'FASTA_Subject':
            return '-subject {}'.format(subject_file)
        elif self.__subject_key == 'DB_BLAST':
            return '-db {}'.format(subject_file)

    def __get_output_name(self, query, subject):
        """
        Generates the default output name.
        :param query: Query file
        :param subject: Subject file
        :return: Output name
        """
        query_basename = os.path.splitext(query.basename)[0]
        subject_basename = os.path.splitext(subject.basename)[0]
        return '{}_{}-{}.{}'.format(self._tool_command, query_basename, subject_basename,
                                    self.__get_output_key().lower())

    def __get_subject_key(self):
        """
        Returns the key of the subject, this can be:
        - FASTA_Subject: FASTA file of the subject sequence
        - DB_BLAST: BLAST database created using makeblastdb.
        :return: Key
        """
        if 'DB_BLAST' and 'FASTA_Subject' in self._tool_inputs:
            raise ValueError("Cannot use DB_BLAST and FASTA_Subject at the same time")
        elif 'DB_BLAST' in self._tool_inputs:
            return 'DB_BLAST'
        elif 'FASTA_Subject' in self._tool_inputs:
            return 'FASTA_Subject'
        else:
            raise ValueError("No subject (FASTA_Subject / DB_BLAST) found")

    def __get_output_key(self):
        """
        Returns the output format key.
        :return: Output key
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

    def _check_command_output(self):
        """
        Checks the command output for errors.
        :return: None
        """
        if 'error' in self._command.stderr.lower():
            raise ValueError("Error executing {}: {}".format(self.name, self._command.stderr.strip()))
