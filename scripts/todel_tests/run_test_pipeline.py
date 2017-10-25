"""
Generates a runtime config yml file for the snakemake pipeline.
Runs the snakemake pipeline.
"""
import os
import yaml

from app.camel import Camel
from app.command.command import Command
from app.io.tooliofile import ToolIOFile
from app.pipeline.snakepipeline import SnakePipeline

DB_LOGGING = False


if __name__ == '__main__':
    # Load the data of the other config file
    with open(os.path.join(os.path.dirname(__file__), 'config.yaml')) as handle:
        config_data = yaml.load(handle)

    # Create a pipeline object
    camel = Camel()
    pipeline = SnakePipeline('GATK somatic calling', camel, DB_LOGGING)

    # Add the job id to the config
    config_data['pipeline_job_id'] = pipeline.job_id
    config_data['pipeline_name'] = pipeline.name
    config_data['logging'] = DB_LOGGING

    # Setting the initial input makes sure that they are logged
    if DB_LOGGING:
        pipeline.set_initial_input({'FASTQ_PE': [ToolIOFile(f) for f in config_data['fastq_pe']]})

    # Create a new config file
    with open('runtime_config.yaml', 'w') as handle:
        yaml.dump(config_data, handle)

    # Execute the snakemake workflow
    command = Command('snakemake --configfile {} --snakefile {}'.format(
        'runtime_config.yaml', os.path.join(os.path.dirname(__file__), 'test.snakefile')
    ))
    command.run_command('.')
    print('Stdout: {}\n'.format(command.stdout))
    print('Stderr: {}\n'.format(command.stderr))
