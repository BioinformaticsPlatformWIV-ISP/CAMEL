import logging
from pathlib import Path
from typing import List, Tuple, Any, Dict, Optional

import yaml

from camel.app.command.command import Command
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.snakemakeexecutionerror import SnakemakeExecutionError
from camel.app.io.tooliofile import ToolIOFile
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
    def add_report_content(report: HtmlReport, report_structure: List[Tuple[str, str, List[str]]]) -> None:
        """
        Adds the content to the HTML report.
        :param report: Report
        :param report_structure: Report structure
        :return: None
        """
        # Add the overview section
        report.add_module_header('Sections')
        section = HtmlReportSection(None)
        overview_list = HtmlElement('ul')
        for title, key, _ in report_structure:
            list_item = HtmlElement('li')
            list_item.add_html_object(HtmlElement('a', title, [('href', '#{}'.format(key))]))
            overview_list.add_html_object(list_item)
        section.add_html_object(overview_list)
        report.add_html_object(section)

        # Add the different sections
        for title, key, items in report_structure:
            report.add_module_header(title, key)
            for pickle in items:
                if len(pickle) == 0:
                    continue
                section = SnakemakeUtils.load_object(pickle)[0].value
                report.add_html_object(section)
                section.copy_files(report.output_dir)
        report.save()

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
    def symlink_input_files(output_dir: Path, file_paths: List[str], file_names: List[str], sanitize: bool = False) ->\
            List[Path]:
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
        if not output_dir.exists():
            output_dir.mkdir(parents=True)
        links = []
        for path, name in zip(file_paths, file_names):
            link_path = output_dir / (FileSystemHelper.make_valid(name) if sanitize else name)
            if link_path.is_symlink():
                link_path.unlink()
            link_path.symlink_to(path)
            links.append(link_path)
        return links

    @staticmethod
    def generate_config_file(config_data: Dict[str, Any], output_dir: Path, output_basename: str = 'config.yml') -> str:
        """
        Generates a configuration file for Snakemake in YAML file format.
        :param config_data: Configuration data
        :param output_dir: Output directory
        :param output_basename: Output basename
        :return: Path to config file
        """
        config_path = output_dir / output_basename
        if not output_dir.exists():
            output_dir.mkdir(parents=True)
        with config_path.open('w') as handle:
            yaml.dump(config_data, handle)
        logging.info(f"Configuration file created: {config_path}")
        return str(config_path)

    @staticmethod
    def run_snakemake(snakefile: str, config_path: str, targets: List[Path], working_dir: Path,
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
        if not working_dir.exists():
            working_dir.mkdir(parents=True)
        command = Command('snakemake {} --snakefile {} --configfile {} --cores {}'.format(
            ' '.join(str(x) for x in targets), snakefile, config_path, threads))
        command.run_command(str(working_dir))
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
            header = f"{informs['_name']} - {informs['_tag']}" if '_tag' in informs else informs['_name']
            section.add_header(header, 3)
            command_txt = informs['_command'].replace(working_dir, '$WORKING')
            section.add_html_object(HtmlElement('code', command_txt, [('class', 'command')]))
        return section

    @staticmethod
    def extracts_fq_input(io_dict: str, key_pe: Optional[str] = 'FASTQ_PE', key_se: Optional[str] = 'FASTQ_SE',
                          drop_se: bool = False) -> Dict[str, List[ToolIOFile]]:
        """
        Extracts a specific FASTQ input dictionary from the standardized FASTQ dictionary.
        :param io_dict: Path to input IO file
        :param key_pe: Key for paired end FASTQ files
        :param key_se: Key for single end FASTQ files
        :param drop_se: If True, SE reads are dropped for PE input
        :return: Reformatted dictionary
        """
        io = SnakemakeUtils.load_object(io_dict)
        if 'PE' in io:
            output_dict = {key_pe: io['PE']}
            se_reads = io.get('SE_FWD', []) + io.get('SE_REV', [])
            if len(se_reads) > 0 and not drop_se:
                output_dict[key_se] = se_reads
        else:
            output_dict = {key_se: io['SE']}
        return output_dict
