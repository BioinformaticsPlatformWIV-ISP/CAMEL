from pathlib import Path

from snakemake.io import expand

from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.snpphylogeny.snakefile.samtools_filtering_all import OUTPUT_FILTERING_ALL


rule filter_variants:
    """
    This rule runs the variant filtering workflow on a sample.
    """
    input:
        VCF = lambda wildcards: config['samples'][wildcards.sample]['VCF'],
        BAM = lambda wildcards: config['samples'][wildcards.sample]['BAM']
    output:
        VCF = Path(config['working_dir'], '{sample}', 'variant_filtering_out.io')
    threads: 4
    params:
        working_dir = lambda wildcards: Path(config['working_dir'], wildcards.sample),
        filtering_options=config['options'],
        sample_name=lambda wildcards: wildcards.sample
    run:
        from camel.app.components.workflows.variantfilteringwrapper import VariantFilteringWrapper
        wrapper = VariantFilteringWrapper(Path(str(params.working_dir), 'variant_filtering'))
        wrapper.run_workflow(str(params.sample_name), Path(str(input.VCF)), Path(str(input.BAM)),
            params.filtering_options, threads)
        SnakemakeUtils.dump_object(wrapper.output, Path(output.VCF))

rule combine_filtering_output:
    """
    Combines the variant filtering output for all samples into a dictionary.
    """
    input:
        VCF = expand(rules.filter_variants.output.VCF, sample=sorted(config['samples'].keys()))
    output:
        IO = Path(config['working_dir'], OUTPUT_FILTERING_ALL)
    params:
        samples = sorted(config['samples'].keys())
    run:
        output_data = {}
        for sample, vcf_file in zip(params.samples, input.VCF):
            output_data[sample] = SnakemakeUtils.load_object(Path(vcf_file))
        SnakemakeUtils.dump_object(output_data, Path(output.IO))
