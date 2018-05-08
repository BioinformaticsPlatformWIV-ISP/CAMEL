import argparse
import logging

import os
import yaml

from camel.app.command.command import Command
from camel.scripts.vtecpipeline import SNAKEFILE_VTEC_MAIN

"""
This test script can be used to run the read trimming component of the VTEC pipeline.
Make sure to set the root directory of your CAMEL installation in the $CAMELPATH environment variable.
Example usage:
run_read_trimming.py
--fastq-pe 4009_1.fastq 4009_2.fastq 
--working-dir /scratch/bebog/working 
--html-out /scratch/bebog/pipeline_out/report.html 
--threads 20
"""


def _parse_arguments():
    """
    Parses the command line arguments.
    :return: Arguments
    """
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--fastq-pe", required=True, nargs=2)
    argument_parser.add_argument("--working-dir", required=True)
    argument_parser.add_argument("--html-out", required=True)
    argument_parser.add_argument("--threads", default=12)
    argument_parser.add_argument("--detection-method", default='fast')
    return argument_parser.parse_args()

if __name__ == '__main__':
    args = _parse_arguments()
    logging.basicConfig(level=logging.DEBUG)

    config_path = os.path.join(args.working_dir, 'config.yml')
    with open(config_path, 'w') as handle:
        yaml.dump({'detection_method': args.detection_method}, handle, default_flow_style=False)
        yaml.dump({'assembler': 'SPAdes'}, handle, default_flow_style=False)
        yaml.dump({'sample_name': 'test_sample'}, handle, default_flow_style=False)
        yaml.dump({'fastq_pe': args.fastq_pe}, handle, default_flow_style=False)
        yaml.dump({'report': args.html_out}, handle, default_flow_style=False)
        yaml.dump({'output_dir': os.path.dirname(args.html_out)}, handle, default_flow_style=False)
        yaml.dump({'working_dir': args.working_dir}, handle, default_flow_style=False)
        yaml.dump({'logging': False}, handle, default_flow_style=False)
        yaml.dump({'pipeline_name': 'Read Trimming'}, handle, default_flow_style=False)
        yaml.dump({'pipeline_job_id': 3}, handle, default_flow_style=False)
        yaml.dump({'skip_assembly': True}, handle, default_flow_style=False)

    command = Command(' '.join([
        '. /home/bebog/venvs/camel_env_3/bin/activate;',
        'export PYTHONPATH=$CAMELPATH;',
        'snakemake',
        f'--cores {args.threads}',
        f'--configfile {config_path}',
        f'--snakefile {SNAKEFILE_VTEC_MAIN}',
    ]))
    command.run_command(args.working_dir)
    print(command.stdout)
    print(command.stderr)
