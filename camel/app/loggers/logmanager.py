import logging
import logging.config
import os

import yaml


class LogManager(object):
    """
    Manages the various logs.
    """
    _step_handlers = []
    _pipeline_handlers = []

    @staticmethod
    def initialize(config_file):
        """
        Initializes the log manager.
        :param config_file: Configuration file
        :return: None
        """
        with open(config_file, 'rt') as f:
            config = yaml.safe_load(f.read())
        logging.config.dictConfig(config)
        LogManager._step_handlers = [h for h in logging.getLogger().handlers if h.get_name().startswith('step')]
        LogManager._pipeline_handlers = [h for h in logging.getLogger().handlers if h.get_name().startswith('pipeline')]
        LogManager.detach_step_handlers()
        LogManager.detach_pipeline_handlers()
        logging.info("Log manager initialized")

    @staticmethod
    def attach_step_handlers(folder):
        """
        Attaches the step handlers and updates the log location.
        :param folder: Folder to store the logs
        :return: None
        """
        for step_handler in LogManager._step_handlers:
            filename = os.path.basename(step_handler.baseFilename)
            step_handler.close()
            step_handler.baseFilename = os.path.join(folder, filename)
            logging.getLogger().addHandler(step_handler)

    @staticmethod
    def detach_step_handlers():
        """
        Detaches the step handlers.
        :return: None
        """
        for step_handler in LogManager._step_handlers:
            if step_handler in logging.getLogger().handlers:
                logging.getLogger().handlers.remove(step_handler)

    @staticmethod
    def attach_pipeline_handlers(folder):
        """
        Attaches the pipeline handlers and updates the log location.
        :param folder: Folder to store the logs
        :return: None
        """
        for pipeline_handler in LogManager._pipeline_handlers:
            filename = os.path.basename(pipeline_handler.baseFilename)
            pipeline_handler.close()
            pipeline_handler.baseFilename = os.path.join(folder, filename)
            logging.getLogger().addHandler(pipeline_handler)

    @staticmethod
    def detach_pipeline_handlers():
        """
        Detaches the pipeline handlers.
        :return: None
        """
        for pipeline_handler in LogManager._pipeline_handlers:
            if pipeline_handler in logging.getLogger().handlers:
                logging.getLogger().handlers.remove(pipeline_handler)
