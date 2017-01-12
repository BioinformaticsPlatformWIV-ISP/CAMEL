import hashlib
import os

import pickle


class FileUtils(object):
    """
    Helper object to perform common operations on the file system.
    """

    @staticmethod
    def hash_file(file_name, block_size=65536):
        """
        Creates a hash for the file with a default block size of 65536 and the sha256 algorithm.
        :param file_name: File that needs to be hashed
        :param block_size: Block size to be used
        :return: String of the hash with alphanumeric symbols
        """
        if not os.path.isfile(file_name):
            raise IOError("'{}' is not a file".format(file_name))
        hasher = hashlib.sha256()
        file_to_hash = open(file_name, 'rb')
        buf = file_to_hash.read(block_size)
        while len(buf) > 0:
            hasher.update(buf)
            buf = file_to_hash.read(block_size)
        return hasher.hexdigest()

    @staticmethod
    def get_all_files(directory_path):
        """
        Returns all files in a directory recursively.
        """
        files_list = []
        for root, directories, files in os.walk(directory_path):
            for file_ in files:
                files_list.append(os.path.join(root, file_))
        return files_list

    @staticmethod
    def hash_directory(path):
        """
        Creates a hash for a folder with a default block size of 65536 and the sha256 algorithm.
        :param path: Directory path
        :return: String of the hash with alphanumeric symbols
        """
        hasher = hashlib.sha256()
        for file_ in sorted(FileUtils.get_all_files(path)):
            hasher.update(FileUtils.hash_file(file_))
        return hasher.hexdigest()

    @staticmethod
    def hash_value(value):
        """
        Creates a hash for a value.
        :param value: Value
        :return: String of the hash with alphanumeric symbols
        """
        hasher = hashlib.sha256()
        hasher.update(pickle.dumps(value))
        return hasher.hexdigest()
