from typing import List

import os
import re


class FileSystemHelper(object):

    @staticmethod
    def make_valid(value: str) -> str:
        """
        Converts arbitrary strings to URL- and filename friendly values.
        :param value: Input value
        :return: URL- and filename friendly value
        """
        value = value.replace(' ', '_')
        return "".join([c for c in value if re.match(r'[\w\-_]', c)])

    @staticmethod
    def get_file_with_extension(input_folder: str, extension: str) -> str:
        """
        Returns a single file with the given extension from the given directory.
        :param input_folder: Input directory
        :param extension: File extension
        :return: Path to file
        """
        all_files = FileSystemHelper.get_files_with_extension(input_folder, extension)
        if len(all_files) == 0:
            raise IOError("No {} file found in '{}'".format(extension, input_folder))
        elif len(all_files) > 1:
            raise IOError("Multiple {} files found in '{}'".format(extension, input_folder))
        return all_files[0]

    @staticmethod
    def get_files_with_extension(folder: str, extension: str) -> List[str]:
        """
        Returns the files with the given extension from the folder.
        :param folder: Input folder
        :param extension: File extension
        :return: List of paths to files
        """
        return [os.path.join(folder, file_) for file_ in os.listdir(folder) if file_.endswith(extension)]
