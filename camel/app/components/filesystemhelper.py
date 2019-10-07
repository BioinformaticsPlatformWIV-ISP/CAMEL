import binascii
import datetime
import logging
from typing import List

import os
import re

from camel.app.command.command import Command


class FileSystemHelper(object):
    """
    This class contains utility function to work with the file system.
    """

    TIMESTAMP_FILENAME = "%Y%m%d-%H%M%S"

    @staticmethod
    def make_valid(value: str) -> str:
        """
        Converts arbitrary strings to URL- and filename friendly values.
        :param value: Input value
        :return: URL- and filename friendly value
        """
        value = value.replace(' ', '_')
        return "".join([c for c in value if re.match(r'[\w\-_\\.]', c)])

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

    @staticmethod
    def get_timestamp_str(timestamp: datetime.datetime = datetime.datetime.now()) -> str:
        """
        Returns the given time stamp as a string that can be used in a filename.
        :param timestamp: Timestamp (default to current time)
        :return: Timestamp as string
        """
        return timestamp.strftime(FileSystemHelper.TIMESTAMP_FILENAME)

    @staticmethod
    def is_gzipped(path: str) -> bool:
        """
        Checks if the given file is compressed with gzip.
        :param path: Path
        :return: True if gzipped, False otherwise
        """
        with open(path, 'rb') as handle:
            magic_number = binascii.hexlify(handle.read(2))
        return magic_number == b'1f8b'

    @staticmethod
    def gzip_extract(input_gz_file: str, output_gz_file) -> None:
        """
        Extracts a GZIP compressed file, the original file is left untouched.
        :param input_gz_file: Input GZ file
        :param output_gz_file: Output path
        :return: None
        """
        logging.info(f"Extracting: {input_gz_file}")
        command = Command(f'gunzip -k -c {input_gz_file} > {output_gz_file}')
        command.run_command('.')
        if not command.returncode == 0:
            raise RuntimeError(f"Cannot extract '{input_gz_file}': {command.stderr}")
