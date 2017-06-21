import os
import unicodedata

import re


class FileSystemHelper(object):
    @staticmethod
    def make_valid(value):
        """
        Converts arbitrary strings to URL- and filename friendly values.
        :param value: Input value
        :return: URL- and filename friendly value
        """
        value = unicodedata.normalize('NFKD', str(value)).encode('ascii', 'ignore')
        value = str(re.sub('[^\w\s-]', '', value).strip())
        value = str(re.sub('[-\s]+', '_', value))
        return str(value)

    @staticmethod
    def get_file_with_extension(input_folder, extension):
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
    def get_files_with_extension(folder, extension):
        """
        Returns the files with the given extension from the folder.
        :param folder: Input folder
        :param extension: File extension
        :return: List of paths to files
        """
        return [os.path.join(folder, file_) for file_ in os.listdir(folder) if file_.endswith(extension)]
