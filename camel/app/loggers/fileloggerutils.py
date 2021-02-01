import gzip
import logging
from pathlib import Path
from typing import Optional

import shutil
import socket

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper


def store_config_file(
        config_file: Path, basename: str, galaxy_job_id: Optional[str] = None, dir_: Path = None) -> Path:
    """
    Exports the config file from the pipeline.
    :param config_file: Config file to export
    :param basename: Basename for tool / pipeline
    :param galaxy_job_id: (Optional) Galaxy job id
    :param dir_: (Optional) Directory to store file, defaults to value from config
    :return: Path to config file
    """
    # Determine output directory
    output_dir = dir_ if dir_ is not None else Path(Camel.get_instance().config.get('logging', {}).get(
        'dir_snakemake_configs'))
    if not output_dir.exists():
        raise RuntimeError(f'Logging directory does not exist: {output_dir}')

    # Determine output path
    export_path = output_dir / '{}.yml.gz'.format('__'.join([
        'config',
        FileSystemHelper.make_valid(basename),
        FileSystemHelper.make_valid(socket.gethostname()),
        FileSystemHelper.get_timestamp_str(),
        galaxy_job_id if galaxy_job_id is not None else 'NA',
    ]))

    # Save config file
    with gzip.open(export_path, 'wb') as file_out, config_file.open('rb') as file_in:
        shutil.copyfileobj(file_in, file_out)
    logging.info(f"Config file stored in: '{export_path}'")
    return export_path


def store_log_file(log_file: Path, basename: str, galaxy_job_id: Optional[str] = None, is_error_log: bool = False,
                   dir_: Optional[Path] = None) -> Path:
    """
    Stores a log file on disk.
    :param log_file: Log file to store
    :param basename: Basename
    :param galaxy_job_id: (Optional) Galaxy job id
    :param is_error_log: Boolean to indicate if this is an error log
    :param dir_: (Optional) Directory to store file, defaults to value from config
    :return: Path to log file
    """
    # Determine output directory
    key_dir = 'dir_error_logs' if is_error_log else 'dir_camel_logs'
    output_dir = dir_ if dir_ is not None else Path(Camel.get_instance().config.get('logging', {}).get(key_dir))
    if not output_dir.exists():
        raise RuntimeError(f'Logging directory does not exist: {output_dir}')

    # Determine output file
    prefix = 'error' if is_error_log else 'camel'
    output_path = output_dir / '{}.txt.gz'.format('__'.join([
        prefix,
        FileSystemHelper.make_valid(basename).lower(),
        FileSystemHelper.make_valid(socket.gethostname()),
        FileSystemHelper.get_timestamp_str(),
        galaxy_job_id if galaxy_job_id is not None else 'NA'
    ]))

    # Save log file
    with gzip.open(output_path, 'wb') as file_out, log_file.open('rb') as file_in:
        shutil.copyfileobj(file_in, file_out)
    logging.debug(f"Log file stored in: '{output_path}'")
    return output_path
