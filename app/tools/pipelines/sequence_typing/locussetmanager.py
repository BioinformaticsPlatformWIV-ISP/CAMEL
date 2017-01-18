import json
import os

import logging

from app.components.filesystemhelper import FileSystemHelper
from app.io.tooliodirectory import ToolIODirectory
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class LocusSetManager(Tool):
    """
    Tool that manages sequence typing locus sets.
    The input is a locus set directory, it is split in a list of nucleotide databases and protein databases. Also
    reports locus set metadata information in the informs.

    The input directory should have the following structure:
    > Scheme
      scheme_metadata.txt
      > Locus_A
        locus_A.fasta
        locus_metadata.txt
      > Locus_B
        locus_B.fasta
        locus_metadata.txt

    Informs:
    - Each FASTA file is a key in the informs, with as value a dictionary containing:
        URL, type, name, allele_id_regex.
    - Locus set metadata, a dictionary containing:
        last updated, description, name.
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :return: None
        """
        super(LocusSetManager, self).__init__('Sequence Typing: Locus Set Manager', '0.1', camel)
        self._tool_outputs['DB_Nucl'] = []
        self._tool_outputs['DB_Prot'] = []

    def _execute_tool(self):
        """
        Run this tool.
        :return: None
        """
        self._informs['gene_names'] = []
        for locus_folder in self.__get_locus_folders():
            fasta_file = FileSystemHelper.get_file_with_extension(locus_folder, 'fasta')
            locus_metadata = LocusSetManager.__get_locus_metadata(locus_folder)
            self._informs[fasta_file] = locus_metadata
            self._informs['gene_names'].append(locus_metadata['name'])
            self.__add_locus_to_output(fasta_file, locus_metadata['type'])
        self.__add_locus_set_metadata()

        profile_definitions = self.__get_sequence_type_profiles()
        if profile_definitions:
            self._informs['has_profile_definitions'] = True
            logging.info("Sequence type definitions found")
            self._tool_outputs['TSV'] = [ToolIOFile(profile_definitions)]
        else:
            self._informs['has_profile_definitions'] = False
        logging.info('{} nucleotide loci found'.format(len(self._tool_outputs['DB_Nucl'])))
        logging.info('{} protein loci found'.format(len(self._tool_outputs['DB_Prot'])))

    def _check_input(self):
        """
        Checks whether the input is correct.
        :return: None
        """
        if 'DIR' not in self._tool_inputs:
            raise ValueError("No 'DIR' input found.")
        if not isinstance(self._tool_inputs['DIR'][0], ToolIODirectory):
            raise IOError("'{}' is not a directory".format(self._tool_inputs['DIR'][0].path))
        super(LocusSetManager, self)._check_input()

    def __add_locus_set_metadata(self):
        """
        Adds the locus set metadata to the informs.
        :return: None
        """
        folder = self._tool_inputs['DIR'][0].path
        try:
            metadata_file = os.path.join(folder, 'scheme_metadata.txt')
            with open(metadata_file) as metadata_handle:
                self._informs['locus_set_metadata'] = json.load(metadata_handle)
        except IOError as err:
            logging.warning('Problem retrieving metadata for {} ({})'.format(folder, err.message))

    def __get_locus_folders(self):
        """
        Returns all the locus folders of this locus set.
        :return: List of paths
        """
        locus_folders = []
        locus_set_folder = self._tool_inputs['DIR'][0]
        logging.info('Checking directory: {}'.format(locus_set_folder))
        for subfolder in os.listdir(locus_set_folder.path):
            locus_folder = os.path.join(locus_set_folder.path, subfolder)
            if os.path.isdir(locus_folder):
                locus_folders.append(os.path.join(locus_set_folder.path, locus_folder))
        if 'limit_locus_sets' in self._parameters:
            return sorted(locus_folders)[0:int(self._parameters['limit_locus_sets'][1])]
        else:
            return sorted(locus_folders)

    @staticmethod
    def __get_locus_metadata(locus_folder):
        """
        Returns the metadata from the given locus folder.
        :return: Dictionary of metadata
        """
        try:
            with open(os.path.join(locus_folder, 'locus_metadata.txt')) as metadata_handle:
                metadata = json.load(metadata_handle)
        except IOError:
            raise IOError("No locus metadata found in '{}'".format(locus_folder))

        locus_informs = {}
        for key in ['name', 'type', 'allele_id_regex', 'url']:
            if key not in metadata:
                raise ValueError("Key '{}' missing from locus metadata".format(key))
            locus_informs[key] = metadata[key]
        return locus_informs

    def __get_sequence_type_profiles(self):
        """
        Returns the sequence type profile definitions.
        :return: Profile definitions file
        """
        try:
            return FileSystemHelper.get_file_with_extension(self._tool_inputs['DIR'][0].path, 'profiles.tsv')
        except IOError:
            return None

    def __add_locus_to_output(self, fasta_file, type_):
        """
        Adds a locus to the output
        :param fasta_file: FASTA file of the locus
        :param type_: Type of the locus ('DNA', 'Peptide')
        :return: None
        """
        if type_ == 'DNA':
            self._tool_outputs['DB_Nucl'].append(ToolIOFile(fasta_file))
        elif type_ == 'peptide':
            self._tool_outputs['DB_Prot'].append(ToolIOFile(fasta_file))
        else:
            raise StandardError("Unrecognized type: '{}' in folder '{}'".format(type_, fasta_file))
