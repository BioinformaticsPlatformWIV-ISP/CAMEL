import logging
from typing import List, Tuple, Any, Dict

import os
import yaml

from camel.app.command.command import Command
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.snakemakeexecutionerror import SnakemakeExecutionError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.snakemake.snakemakeutils import SnakemakeUtils


class SnakePipelineUtils(object):
    """
    This class contains utility functions for Snakemake pipelines.
    """

    DATE_FORMAT = '%d/%m/%Y - %X'

    @staticmethod
    def create_input_section(sample_name: str, date: str, pipeline_version: str, input_files: str,
                             extra_data: List[Tuple[str, str]]) -> HtmlReportSection:
        """
        Creates the input section for the HTML report.
        :param sample_name: Sample name
        :param date: Analysis date
        :param pipeline_version: Pipeline version
        :param input_files: Input files
        :param extra_data: Extra data to include in the input section
        :return:
        """
        table_data = [
            ['Sample:', sample_name],
            ['Analysis date', date],
            ['Pipeline version:', pipeline_version],
            ['Input files:', input_files],
        ]
        for key, value in extra_data:
            table_data.append([f'{key}:', value])
        section = HtmlReportSection('Input')
        section.add_table(table_data, table_attributes=[('class', 'information')])
        return section

    @staticmethod
    def create_sections_overview(report_sections: List[Tuple[str, str, List[str]]]) -> HtmlReportSection:
        """
        Creates the sections overview for the HTML report.
        :param report_sections: Overview of the report sections (name, abbreviation, path)
        :return: HTML report section
        """
        section = HtmlReportSection(None)
        overview_list = HtmlElement('ul')
        for title, key, _ in report_sections:
            list_item = HtmlElement('li')
            list_item.add_html_object(HtmlElement('a', title, [('href', '#{}'.format(key))]))
            overview_list.add_html_object(list_item)
        section.add_html_object(overview_list)
        return section

    @staticmethod
    def create_empty_report_section(title: str, output_file: str, header_level: int = 3) -> None:
        """
        Creates an empty report section.
        :param title: Section title
        :param output_file: Output file
        :param header_level: Header level
        :return: None
        """
        section = HtmlReportSection(title, header_level)
        section.add_paragraph('Analysis disabled')
        SnakemakeUtils.dump_object([ToolIOValue(section)], output_file)

    @staticmethod
    def symlink_input_files(output_dir: str, file_paths: List[str], file_names: List[str], sanitize: bool = False) ->\
            List[str]:
        """
        Creates symlinks with the given names for the given files.
        This can be used for files that come from Galaxy that have a fixed name (dataset_XXXXX.dat).
        :param output_dir: Directory to save symlinks
        :param file_paths: Input file paths
        :param file_names: Input file names
        :param sanitize: If True, file names are sanitized
        :return: List of absolute paths to symlinks
        """
        if len(file_names) != len(file_paths):
            raise ValueError("File names ({}) and file paths ({}) should be the same length".format(
                len(file_names), len(file_paths)))
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        links = []
        for path, name in zip(file_paths, file_names):
            link_path = os.path.join(output_dir, FileSystemHelper.make_valid(name) if sanitize else name)
            if os.path.islink(link_path):
                os.remove(link_path)
            os.symlink(path, link_path)
            links.append(link_path)
        return links

    @staticmethod
    def generate_config_file(config_data: Dict[str, Any], output_dir: str, output_basename: str = 'config.yml') -> str:
        """
        Generates a configuration file for Snakemake in YAML file format.
        :param config_data: Configuration data
        :param output_dir: Output directory
        :param output_basename: Output basename
        :return: Path to config file
        """
        config_path = os.path.join(output_dir, output_basename)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        with open(config_path, 'w') as handle:
            yaml.dump(config_data, handle)
        logging.info(f"Configuration file created: {config_path}")
        return config_path

    @staticmethod
    def run_snakemake(snakefile: str, config_path: str, targets: List[str], working_dir: str,
                      threads: int = 8) -> None:
        """
        Helper function to run snakemake workflows.
        :param snakefile: Workflow snakefile
        :param config_path: Path to configuration file
        :param targets: Target output files
        :param working_dir: Working directory
        :param threads: Number of threads to use
        :return: None
        """
        if not os.path.isdir(working_dir):
            os.makedirs(working_dir)
        command = Command('snakemake --snakefile {} --configfile {} {} --cores {}'.format(
            snakefile, config_path, ' '.join(targets), threads))
        command.run_command(working_dir)
        print(f'- Stdout: -\n{command.stdout}')
        print(f'- Stderr: -\n{command.stderr}')
        if command.returncode != 0:
            raise SnakemakeExecutionError(command.stdout, command.stderr)

    @staticmethod
    def create_commands_section(tool_informs: List[Dict[str, Any]], working_dir: str) -> HtmlReportSection:
        """
        Creates a section with an overview of the commands.
        :param tool_informs: Tool informs
        :param working_dir: Working directory
        :return: Commands section
        """
        section = HtmlReportSection('Commands')
        logging.debug(f"Exporting command for {len(tool_informs)} tools")
        for informs in tool_informs:
            section.add_header(informs['_name'], 3)
            command_txt = informs['_command'].replace(working_dir, '$WORKING')
            section.add_html_object(HtmlElement('code', command_txt, [('class', 'command')]))
        return section
