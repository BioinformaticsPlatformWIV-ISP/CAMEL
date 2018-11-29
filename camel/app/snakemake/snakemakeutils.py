import bz2
import logging
import pickle
import time
from typing import Any, Optional, List

import os
from shutil import copyfileobj

from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class SnakemakeUtils(object):

    """
    This class contains utility functions for working with snakemake and CAMEL.
    """

    @staticmethod
    def dump_object(obj: Any, path: str) -> None:
        """
        Dumps an object in a pickle.
        :param obj: Object to dump
        :param path: Path to store the pickle
        :return: None
        """
        logging.debug("Dumping object '{!r}' in file '{}'".format(obj, path))
        with open(path, 'wb') as handle:
            pickle.dump(obj, handle)

    @staticmethod
    def load_object(path: str) -> Any:
        """
        Loads the object from the given pickle.
        :param path: Path
        :return: Object
        """
        logging.debug("Loading object from file '{}'".format(path))
        with open(path, 'rb') as handle:
            obj = pickle.load(handle)
        logging.debug("'{!r}' loaded".format(obj))
        return obj

    @staticmethod
    def get_io_object(value: Any) -> Any:
        """
        Returns the value as a CAMEL IO object.
        - If it is a file path a ToolIOFile object is returned
        - If it is a directory path a ToolIOFolder object is returned
        - Else a ToolIOValue is returned
        :param value: Input value
        :return: ToolIO object
        """
        if os.path.isfile(value):
            converted_value = ToolIOFile(value)
        elif os.path.isdir(value):
            converted_value = ToolIODirectory(value)
        else:
            converted_value = ToolIOValue(value)
        logging.info("'{}' converted to {!r}".format(value, converted_value))
        return converted_value

    @staticmethod
    def add_pickle_input(tool: Tool, key: str, path: str, optional: bool=False) -> None:
        """
        Adds a pickled input to a tool. For optional input whose value is empty, it is skipped.
        :param tool: Tool
        :param key: Key
        :param path: Pickle path
        :param optional: True for optional input, False otherwise
        :return: None
        """
        logging.debug("Adding pickled input with key '{}' from file '{}' to tool '{}'".format(
            key, path, tool.name))
        value = SnakemakeUtils.load_object(path)
        if optional and len(value) == 0:
            logging.debug("Optional Input '{!r}' empty, skipped".format(key))
        else:
            tool.add_input_files({key: value})

    @staticmethod
    def dump_tool_output(tool: Tool, key: str, path: str) -> None:
        """
        Dumps a tool output to a Camel IO pickle.
        :param tool: Tool
        :param key: Key
        :param path: Pickle path
        :return: None
        """
        logging.debug("Dumping output with key '{}' from tool '{}' to Camel IO pickle '{}'".format(
            key, tool.name, path))
        if key not in tool.tool_outputs:
            raise KeyError("Tool '{}' has no output '{}'".format(tool.name, key))
        SnakemakeUtils.dump_object(tool.tool_outputs[key], path)

    @staticmethod
    def add_pickle_inputs(tool: Tool, snake_input: Any, keys: Optional[List[str]]=None,
                          excluded_keys: Optional[List[str]]=None, optionals: Optional[List[str]]=None) -> None:
        """
        Adds pickled inputs from the snakemake input. If 'optionals' is specified, any optional input in that
        list will be skipped if its value is empty (no input file).
        :param tool: Tool
        :param snake_input: Snakemake input
        :param keys: Keys to add. If None, all keys are added
        :param excluded_keys: For the keys in this list the files are not added
        :param optionals: list of keys specifying optional inputs
        :return: None
        """
        logging.info("Adding pickled inputs from snakemake input")
        if optionals is None:
            optionals = []
        if keys is None:
            keys = snake_input.keys()
        for key in keys:
            logging.debug("Adding input '{}'".format(key))
            if key not in snake_input.keys():
                raise KeyError("Key '{}' not found in snakemake input".format(key))
            if (excluded_keys is not None) and (key in excluded_keys):
                continue
            with open(snake_input[key], 'rb') as handle:
                value = pickle.load(handle)
            if key.startswith('INFORMS'):
                inform_key = '_'.join(key.split('_')[1:])
                tool.add_input_informs({inform_key: value})
                logging.debug("Informs '{!r}' added".format(value))
            else:
                if key in optionals and len(value) == 0:
                    logging.debug("Optional Input '{!r}' empty, skipped".format(key))
                    continue
                tool.add_input_files({key: value})
                logging.debug("Input '{!r}' added".format(value))

    @staticmethod
    def dump_tool_outputs(tool: Tool, snake_output: Any, keys: Optional[List[str]]=None,
                          ignore_missing_output=False) -> None:
        """
        Dumps the tool outputs in pickles.
        :param tool: Tool
        :param snake_output: Snake output
        :param keys: Keys to dump
        :param ignore_missing_output: If False, an error is raised when an output is not generated
        :return: None
        """
        logging.info("Dumping tool outputs")
        if keys is None:
            keys = snake_output.keys()
        for key in keys:
            if key in tool.tool_outputs:
                with open(snake_output[key], 'wb') as handle:
                    pickle.dump(tool.tool_outputs[key], handle)
            elif key == 'INFORMS':
                with open(snake_output[key], 'wb') as handle:
                    pickle.dump(tool.informs, handle)
            else:
                message = "Output '{}' not generated".format(key)
                if ignore_missing_output is True:
                    logging.warning(message)
                else:
                    raise ValueError(message)

    @staticmethod
    def pickle_snake_input(snake_input: Any, snake_output: Any, keys: Optional[List[str]]=None) -> None:
        """
        Converts snakemake input to CAMEL IO pickles.
        For every key, it will attempt to convert every value to an IO object (see IO object function) and store the
        generated pickle in the corresponding file specified in the snake output.
        :param snake_input: Snake input
        :param snake_output: Snake output
        :param keys: If specified, only those keys are converted.
        :return: None
        """
        logging.info("Converting snake input '{!r}' to pickles".format(snake_input))
        if keys is None:
            keys = snake_input.keys()
        for key in keys:
            if key not in snake_input.keys():
                raise KeyError("Key '{}' not found in snakemake input".format(key))
            input_list = snake_input[key]
            if key not in snake_output.keys():
                raise ValueError("Output key '{}' not found.".format(key))
            list_io_objects = [SnakemakeUtils.get_io_object(i) for i in input_list]
            SnakemakeUtils.dump_object(list_io_objects, snake_output[key])

    @staticmethod
    def run_tool(tool: Tool, snake_input: Any, snake_output: Any, working_dir: str) -> None:
        """
        Runs a tool and collects / converts the output and input.
        :param tool: Tool
        :param snake_input: Snakemake input
        :param snake_output: Snakemake output
        :param working_dir: Working directory
        :return: None
        """
        logging.info("Running tool: {}".format(tool.name))
        SnakemakeUtils.add_pickle_inputs(tool, snake_input)
        tool.run(working_dir)
        SnakemakeUtils.dump_tool_outputs(tool, snake_output)

    @staticmethod
    def dump_config(config_path: str, job_id: int, folder: str) -> None:
        """
        Dumps the Snakemake config file in the given folder.
        :param config_path: Path of the config file
        :param job_id: Pipeline job id
        :param folder: Folder to store the config file
        :return: None
        """
        logging.info("Dumping config file '{}' in: {}".format(os.path.basename(config_path), folder))
        outfile = os.path.join(folder, 'config_{}_{}.bz2'.format(job_id, time.strftime("%Y%m%d")))
        with open(config_path, 'rb') as handle_in:
            with bz2.BZ2File(outfile, 'wb', compresslevel=9) as handle_out:
                copyfileobj(handle_in, handle_out)
