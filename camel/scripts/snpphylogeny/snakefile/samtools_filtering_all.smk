import os

from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.snpphylogeny.snakefile.samtools_filtering_all import OUTPUT_FILTERING_ALL


rule Run_variant_filtering_workflow:
    """
    This rule runs the variant filtering workflow on a sample.
    """
    input:
        VCF=lambda wildcards: config['samples'][wildcards.sample]['VCF'],
        BAM=lambda wildcards: config['samples'][wildcards.sample]['BAM']
    output:
        VF_OUTPUT=os.path.join(config['working_dir'], '{sample}', 'variant_filtering_out.io')
    threads: 4
    params:
        working_dir=os.path.join(config['working_dir'], '{sample}'),
        filtering_options=config['options'],
        sample_name=lambda wildcards: wildcards.sample
    run:
        from camel.app.components.workflows.variantfilteringwrapper import VariantFilteringWrapper
        wrapper = VariantFilteringWrapper(os.path.join(params.working_dir, 'variant_filtering'))
        wrapper.run_workflow(input.VCF, input.BAM, params.filtering_options, threads)
        SnakemakeUtils.dump_object(wrapper.output, output.VF_OUTPUT)

rule Collect_variant_filtering_output:
    """
    Combines the variant filtering output for all samples into a dictionary.
    """
    input:
        VF_OUTPUT=expand(os.path.join(config['working_dir'], '{sample}', 'variant_filtering_out.io'),
                         sample=sorted(config['samples'].keys()))
    output:
        os.path.join(config['working_dir'], OUTPUT_FILTERING_ALL)
    params:
        samples=sorted(config['samples'].keys())
    run:
        output_data = {}
        for i in range(0, len(params.samples)):
            output_data[params.samples[i]] = SnakemakeUtils.load_object(input.VF_OUTPUT[i])
        SnakemakeUtils.dump_object(output_data, output[0])
