import logging


class LogManager(object):
    """
    Manages the various logs.
    """

    _formatter = logging.Formatter('CAMEL - %(levelname)s - %(asctime)s - %(message)s')

    @staticmethod
    def get_file_handler(location, level):
        """
        Returns a file handler.
        :param location: Log file location
        :param level: logging level
        :return: File handler
        """
        file_handler = logging.FileHandler(location)
        file_handler.setLevel(level)
        file_handler.setFormatter(LogManager._formatter)
        return file_handler
