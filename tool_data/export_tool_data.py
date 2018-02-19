import argparse

import logging
import os
import yaml

from app.camel import Camel
from app.services.dbtoolservice import DbToolService


"""
This scripts is used to export tool data from the database to a YAML file.
The YAML file name is generated based on the name and version of the tool. 

It can be ran from the command line:
export_tool_data.py --tool-name TOOL_NAME --version VERSION

It can also be ran from within Python by importing it:
from tool_data import export_tool_data
export_tool_data(camel, 'name', 'version')
"""


def _parse_arguments():
    """
    Parses the command line arguments.
    :return: Arguments
    """
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('--tool-name', required=True)
    argument_parser.add_argument('--tool-version', required=True)
    argument_parser.add_argument('--force', required=False, action='store_true')
    return argument_parser.parse_args()


def _get_tool_data(camel: Camel, tool_name: str, version: str):
    """
    Returns the data for the given tool.
    :param camel: CAMEL instance
    :param tool_name: Tool name
    :param version: Tool version
    :return: Tool data
    """
    tool_service = DbToolService(tool_name, version, camel.connection)

    # Collect tool data
    tool_data = {
        'tool_command': tool_service.get_tool_command(),
        'dependencies': tool_service.get_dependencies()
    }
    mandatory_parameters = tool_service.get_names_mandatory_parameter()
    default_parameters = tool_service.get_default_parameters()

    # Add the other parameters
    tool_data['parameters'] = {}
    for key, parameter in tool_service.get_all_parameters().items():
        tool_data['parameters'][parameter.name] = {
            'option': parameter.option,
            'value': parameter.value,
            'mandatory': parameter.name in mandatory_parameters,
            'default': parameter.name in default_parameters
        }
    return tool_data


def export_tool_data(camel: Camel, tool_name: str, tool_version: str, force=False) -> None:
    """
    Exports the data for a given tool.
    :param camel: CAMEL instance
    :param tool_name: Tool name
    :param tool_version: Tool version
    :param force: If True, existing files are overwritten
    :return: None
    """
    output_path = Camel.get_tool_data_path(tool_name, tool_version)
    logging.info("Exporting tool data for {} {} to: {}".format(tool_name, tool_version, output_path))
    if os.path.exists(output_path) and force is False:
        raise FileExistsError("Tool data file already exists: {}".format(output_path))
    tool_data = _get_tool_data(camel, tool_name, tool_version)
    with open(output_path, 'w') as handle:
        yaml.dump(tool_data, handle)


if __name__ == '__main__':
    c = Camel()
    args = _parse_arguments()
    export_tool_data(c, args.tool_name, args.tool_version, args.force)
