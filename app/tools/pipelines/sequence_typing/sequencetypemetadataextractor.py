import logging
import os
from collections import OrderedDict

from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class SequenceTypeMetadataExtractor(Tool):
    """
    Tool that extracts the relevant metadata from the MLST scheme based on the detected sequence type.
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: Camel instance
        :return: None
        """
        super(SequenceTypeMetadataExtractor, self).__init__('Sequence Typing: Metadata Extractor', '0.1', camel)

    def _execute_tool(self):
        self.__add_st_metadata_to_informs(self._tool_inputs['TSV'][0])
        logging.info("Metadata for sequence type: {}".format(self.informs['metadata'].items()))
        filename = os.path.join(self._folder, 'metadata.tsv')
        self.__save_sequence_type_metadata(filename)
        self._tool_outputs['TSV'] = [ToolIOFile(filename)]

    def _check_input(self):
        """
        Checks whether the input is correct.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise ValueError("No sequence type definitions 'TSV' input found.")
        if 'type_detection' not in self._input_informs:
            raise ValueError("No sequence type found.")
        if 'database_info' not in self._input_informs:
            raise ValueError("No database info found")

    @staticmethod
    def __get_sequence_type_line(scheme_file, sequence_type, index=0):
        """
        Returns the line corresponding to the given sequence type.
        :param scheme_file: Scheme file
        :param sequence_type: Sequence type
        :param index: Index of the sequence type
        :return: Line containing the sequence type
        """
        clean_sequence_type = sequence_type.replace('?', '').replace('*', '')
        with open(scheme_file.path) as input_mlst:
            lines = input_mlst.readlines()
            for line in lines:
                st = line.split('\t')[index]
                if st == clean_sequence_type:
                    return line

    @staticmethod
    def __get_sequence_type_header(scheme_file):
        """
        Returns the header
        :param scheme_file:
        :return:
        """
        with open(scheme_file.path) as input_mlst:
            return input_mlst.readline().strip()

    def __add_st_metadata_to_informs(self, mlst_scheme_file):
        """
        Extracts the metadata from the input files and adds them to the metadata informs.
        :param mlst_scheme_file: File containing the MLST scheme
        :return: None
        """
        self.informs['metadata'] = OrderedDict()
        metadata_columns = self.__get_metadata_columns(mlst_scheme_file)
        sequence_type = self._input_informs['type_detection']['sequence_type']
        scheme_line = self.__get_sequence_type_line(mlst_scheme_file, sequence_type)
        for index, name in sorted(metadata_columns.iteritems()):
            if scheme_line:
                value = scheme_line.split('\t')[index].strip()
                if value == '':
                    self.informs['metadata'][name] = '-'
                else:
                    self.informs['metadata'][name] = value
            else:
                self.informs['metadata'][name] = '-'

    def __get_metadata_columns(self, mlst_scheme_file):
        """
        Returns the metadata columns.
        :return: Metadata columns
        """
        header = self.__get_sequence_type_header(mlst_scheme_file).split('\t')
        gene_names = self._input_informs['database_info']['gene_names']
        metadata_columns = {}
        for i in range(0, len(header)):
            if header[i] not in gene_names:
                metadata_columns[i] = header[i].strip()
        return metadata_columns

    def __save_sequence_type_metadata(self, filename):
        """
        Writes the informs to a file.
        :param filename: Filename
        :return: None
        """
        metadata = self.informs['metadata']
        with open(filename, 'w') as output_file:
            output_file.write('\t'.join(metadata.keys()))
            output_file.write('\n')
            output_file.write('\t'.join(metadata.values()))
            output_file.write('\n')
