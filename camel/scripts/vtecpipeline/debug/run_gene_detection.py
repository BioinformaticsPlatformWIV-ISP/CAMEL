import argparse
import logging

import os
import yaml

from camel.app.command.command import Command
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.vtecpipeline import SNAKEFILE_VTEC_MAIN

"""
This test script can be used to run the gene detection component of the VTEC pipeline. 
Make sure to set the root directory of your CAMEL installation in the $CAMELPATH environment variable.
Example usage (alignment based):
run_gene_detection.py
--fasta ecoli_O104-H4.fasta
--working-dir /scratch/bebog/working
--html-out /scratch/bebog/pipeline_out/report.html
--threads 20
--resistance ARG-ANNOT,ResFinder
--virulence Shiga-Toxin_genes

Example usage (read mapping based):
run_gene_detection.py
--fastq-pe 4009_1.fastq 4009_2.fastq
--detection-method normal
--working-dir /scratch/bebog/working
--html-out /scratch/bebog/pipeline_out/report.html
--threads 20
--resistance ARG-ANNOT,ResFinder
--virulence Shiga-Toxin_genes
"""


def _parse_arguments():
    """
    Parses the command line arguments.
    :return: Arguments
    """
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--fasta")
    argument_parser.add_argument("--fastq-pe", nargs=2, default=[])
    argument_parser.add_argument("--working-dir", required=True)
    argument_parser.add_argument("--resistance")
    argument_parser.add_argument("--virulence")
    argument_parser.add_argument("--serotype")
    argument_parser.add_argument("--plasmid")
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
        yaml.dump({'pipeline_name': 'VTEC Pipeline'}, handle, default_flow_style=False)
        yaml.dump({'pipeline_job_id': 3}, handle, default_flow_style=False)
        db_config = {}
        if args.resistance is not None:
            db_config['resistance'] = args.resistance.split(',')
        if args.virulence is not None:
            db_config['virulence'] = args.virulence.split(',')
        if args.serotype is not None:
            db_config['serotype'] = args.serotype.split(',')
        if args.plasmid is not None:
            db_config['plasmid'] = args.plasmid.split(',')
        yaml.dump({'gene_detection': db_config}, handle, default_flow_style=False)
        yaml.dump({'skip_assembly': True}, handle, default_flow_style=False)
        yaml.dump({'skip_trimming': True}, handle, default_flow_style=False)

    # Create the assembly file
    if args.fasta is not None:
        if not os.path.isdir(os.path.join(args.working_dir, 'assembly')):
            os.makedirs(os.path.join(args.working_dir, 'assembly'))
        SnakemakeUtils.dump_object([ToolIOFile(args.fasta)], os.path.join(args.working_dir, 'assembly', 'fasta.io'))

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
