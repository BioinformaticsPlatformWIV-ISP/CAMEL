import os

from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.snpphylogeny.snakefile.trimming_all import TRIMMING_ALL

rule Run_trimming_workflow:
    """
    This rule runs the trimming workflow on a sample.
    """
    input:
        FASTQ_PE=lambda wildcards: config['samples'][wildcards.sample]
    output:
        TRIMMING_OUTPUT=os.path.join(config['working_dir'], '{sample}', 'trimming', 'trimming_out.io')
    threads: 4
    params:
        working_dir=os.path.join(config['working_dir'], '{sample}')
    run:
        from camel.app.components.workflows.readtrimmingwrapper import ReadTrimmingWrapper
        trimming = ReadTrimmingWrapper(os.path.join(params.working_dir, 'trimming'))
        trimming.run_workflow(input.FASTQ_PE, threads)
        SnakemakeUtils.dump_object(trimming.output, output.TRIMMING_OUTPUT)

rule Collect_trimming_output:
    """
    Combines the trimming output for all samples into a dictionary.
    """
    input:
        TRIMMING_OUTPUT=expand(os.path.join(config['working_dir'], '{sample}', 'trimming', 'trimming_out.io'),
                               sample=sorted(config['samples'].keys()))
    output:
        os.path.join(config['working_dir'], TRIMMING_ALL)
    params:
        samples=sorted(config['samples'].keys())
    run:
        output_data = {}
        for i in range(0, len(params.samples)):
            output_data[params.samples[i]] = SnakemakeUtils.load_object(input.TRIMMING_OUTPUT[i])
        SnakemakeUtils.dump_object(output_data, output[0])
