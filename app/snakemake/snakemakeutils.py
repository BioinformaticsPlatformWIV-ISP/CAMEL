import logging
import pickle

import os

from app.io.tooliodirectory import ToolIODirectory
from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue


def dump_object(obj, path):
    """
    Dumps an object in a pickle.
    :param obj: Object to dump
    :param path: Path to store the pickle
    :return: None
    """
    logging.debug("Dumping object '{!r}' in file '{}'".format(obj, path))
    pickle.dump(obj, open(path, 'wb'))


def load_object(path):
    """
    Loads the object from the given pickle.
    :param path: Path
    :return: Object
    """
    logging.debug("Loading object from file '{}'".format(path))
    obj = pickle.load(open(path, 'rb'))
    logging.debug("'{!r}' loaded".format(obj))
    return obj


def get_io_object(value):
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


def add_pickle_input(tool, key, path):
    """
    Adds a pickled input to a tool.
    :param tool: Tool
    :param key: Key
    :param path: Pickle path
    :return: None
    """
    logging.debug("Adding pickled input with key '{}' from file '{}' to tool '{}'".format(
        key, path, tool.name))
    tool.add_input_files({key: load_object(path)})


def dump_tool_output(tool, key, path):
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
    dump_object(tool.tool_outputs[key], path)


def add_pickle_inputs(tool, snake_input, keys=None):
    """
    Adds pickled inputs from the snakemake input.
    :param tool: Tool
    :param snake_input: Snakemake input
    :param keys: Keys to add. If None, all keys are added
    :return: None
    """
    logging.info("Adding pickled inputs from snakemake input")
    if keys is None:
        keys = snake_input.keys()
    for key in keys:
        logging.debug("Adding input '{}'".format(key))
        if key not in snake_input.keys():
            raise KeyError("Key '{}' not found in snakemake input".format(key))
        value = pickle.load(open(snake_input[key], 'rb'))
        if key.startswith('INFORMS'):
            inform_key = ''.join(key.split('_')[1:])
            tool.add_input_informs({inform_key: value})
            logging.debug("Informs '{!r}' added".format(value))
        else:
            tool.add_input_files({key: value})
            logging.debug("Input '{!r}' added".format(value))


def dump_tool_outputs(tool, snake_output, keys=None, ignore_missing_output=False):
    """
    Dumps the tool outputs in pickles.
    :param tool: Tool
    :param snake_output: Snake output
    :param keys: Keys to dump
    :param ignore_missing_output: If False, an error is raised when an output is not generated
    :return: None
    """
    if keys is None:
        keys = snake_output.keys()
    for key in keys:
        if key in tool.tool_outputs:
            pickle.dump(tool.tool_outputs[key], open(snake_output[key], 'wb'))
        elif key == 'INFORMS':
            pickle.dump(tool.informs, open(snake_output[key], 'wb'))
        else:
            message = "Output '{}' not generated".format(key)
            if ignore_missing_output is True:
                logging.warning(message)
            else:
                raise ValueError(message)


def pickle_snake_input(snake_input, snake_output, keys=None):
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
        list_io_objects = [get_io_object(i) for i in input_list]
        dump_object(list_io_objects, snake_output[key])


def run_tool(tool, snake_input, snake_output, working_dir):
    """
    Runs a tool and collects / converts the output and input.
    :param tool: Tool
    :param snake_input: Snakemake input
    :param snake_output: Snakemake output
    :param working_dir: Working directory
    :return: None
    """
    add_pickle_inputs(tool, snake_input)
    tool.run(working_dir)
    dump_tool_outputs(tool, snake_output)
