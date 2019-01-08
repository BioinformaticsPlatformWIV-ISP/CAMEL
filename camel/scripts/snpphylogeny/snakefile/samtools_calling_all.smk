import os

from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.snpphylogeny.snakefile.samtools_calling_all import OUTPUT_CALLING_ALL

rule Run_variant_calling_workflow:
    """
    This rule runs the variant calling workflow on a sample.
    """
    output:
        VC_OUTPUT=os.path.join(config['working_dir'], '{sample}', 'variant_calling_out.io')
    threads: 4
    params:
        working_dir=os.path.join(config['working_dir'], '{sample}'),
        ref_info=config['reference_info'],
        calling_options=config['options'],
        sample_name=lambda wildcards: wildcards.sample,
        sample_config=lambda wildcards: config['samples'][wildcards.sample]
    run:
        from camel.app.components.workflows.variantcallingwrapper import VariantCallingWrapper
        input_files = VariantCallingWrapper.VariantCallingInput(
            pe_reads=[ToolIOFile(x) for x in params.sample_config['PE']],
            se_reads_fwd=ToolIOFile(params.sample_config['SE_FWD']) if 'SE_FWD' in params.sample_config else None,
            se_reads_rev=ToolIOFile(params.sample_config['SE_REV']) if 'SE_REV' in params.sample_config else None
        )
        wrapper = VariantCallingWrapper(os.path.join(params.working_dir, 'variant_calling'))
        wrapper.run_workflow(params.ref_info, params.sample_name, input_files, params.calling_options, threads)
        SnakemakeUtils.dump_object(wrapper.output, output.VC_OUTPUT)

rule Collect_variant_calling_output:
    """
    Combines the variant calling output for all samples into a dictionary.
    """
    input:
        VC_OUTPUT=expand(os.path.join(config['working_dir'], '{sample}', 'variant_calling_out.io'),
                         sample=sorted(list(config['samples'].keys())))
    output:
        os.path.join(config['working_dir'], OUTPUT_CALLING_ALL)
    params:
        samples=sorted(config['samples'].keys())
    run:
        output_data = {}
        for i in range(0, len(params.samples)):
            output_data[params.samples[i]] = SnakemakeUtils.load_object(input.VC_OUTPUT[i])
        SnakemakeUtils.dump_object(output_data, output[0])
