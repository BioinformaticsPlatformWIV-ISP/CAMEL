"""
This Snakefile performs variant calling using the clair3 pipeline.
"""
from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import variant_calling_clair3

camel = Camel.get_instance()


rule clair3_variant_calling_prep_reference:
    """
    Converts the reference to a Snakemake / CAMEL compatible format. Creates the reference metadata.
    """
    output:
        FASTA = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_FASTA,
        INFORMS = Path(config['working_dir']) / 'variant_calling' / 'reference' / 'informs.io'
    params:
        reference = config['variant_calling'].get('reference')
    run:
        from camel.app.io.tooliofile import ToolIOFile
        SnakemakeUtils.dump_object([ToolIOFile(Path(params.reference['path']))], Path(output.FASTA))
        SnakemakeUtils.dump_object(params.reference, Path(output.INFORMS))

rule clair3_variant_calling_alignment_sorting:
    """
    Sorts the alignment.
    """
    input:
        BAM = Path(config['working_dir']) / variant_calling_clair3.INPUT_VARIANT_CALLING_BAM
    output:
        BAM = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_BAM
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'alignment_sorting'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort(camel)
        step = Step(str(rule), samtools_sort, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_sort, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_sort, output)

rule clair3_variant_calling_index_bam:
    """
    Index the bam file.
    """
    input:
        BAM = rules.clair3_variant_calling_alignment_sorting.output.BAM
    output:
        BAM = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_BAM_INDEX
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'alignment_sorting'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        samtools_index = SamtoolsIndex(camel)
        step = Step(str(rule), samtools_index, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_index, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_index, output)

rule clair3_variant_calling_index_fasta:
    """
    Index the fasta file.
    """
    input:
        FASTA = rules.clair3_variant_calling_prep_reference.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_FASTA_INDEX
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'reference'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools_fastaindex = SamtoolsFastaIndex(camel)
        step = Step(str(rule), samtools_fastaindex, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_fastaindex, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_fastaindex, output)

rule clair3_variant_calling:
    """
    This step actually calls variants with clair3
    """
    input:
        FASTA = rules.clair3_variant_calling_index_fasta.output.FASTA,
        BAM = rules.clair3_variant_calling_alignment_sorting.output.BAM,
        BAM_INDEX = rules.clair3_variant_calling_index_bam.output.BAM
    output:
        VCF = Path(config['working_dir']) / 'variant_calling' / 'calling' / 'vcf_gz.io',
        INFORMS = Path(config['working_dir']) / 'variant_calling' / 'calling' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'calling',
        long_indel = config['variant_calling']['long_indel'],
        haploid_precise = config['variant_calling']['haploid_precise'],
        include_ctgs = config['variant_calling']['include_ctgs'],
        no_phasing = config['variant_calling']['no_phasing'],
        model_path = config['model_path'],
        platform = config['variant_calling']['platform']
    threads: 8
    run:
        from camel.app.tools.clair3.clair3 import Clair3
        clair3 = Clair3(camel)
        SnakemakeUtils.add_pickle_inputs(clair3,input)
        step = Step(str(rule), clair3, camel, params.running_dir, config)
        clair3.update_parameters(model_path=params.model_path, platform=params.platform,
            output_path=str(params.running_dir), threads=threads)
        if params.long_indel:
            clair3.update_parameters(long_indel='')
        if params.haploid_precise:
            clair3.update_parameters(haploid_precise='')
        if params.include_ctgs:
            clair3.update_parameters(include_ctgs='')
        if params.no_phasing:
            clair3.update_parameters(no_phasing='')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(clair3, output)

rule clair3_variant_calling_normalize_indels:
    """
    Normalizes indels.
    """
    input:
        VCF_GZ = rules.clair3_variant_calling.output.VCF,
        FASTA = rules.clair3_variant_calling_prep_reference.output.FASTA
    output:
        VCF_GZ = Path(config['working_dir']) / 'variant_calling' / 'norm' / 'vcf_gz.io'
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'norm'
    run:
        from camel.app.tools.bcftools.bcftoolsnorm import BcftoolsNorm
        bcftools_norm = BcftoolsNorm(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(bcftools_norm, input)
        step = Step(str(rule), bcftools_norm, camel, params.running_dir, config)
        bcftools_norm.update_parameters(output_format='z')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_norm, output)

rule clair3_variant_calling_index_vcf_gz:
    """
    Indexes the VCF file.
    """
    input:
        VCF_GZ = rules.clair3_variant_calling_normalize_indels.output.VCF_GZ
    output:
        VCF_GZ = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'norm'
    run:
        from camel.app.tools.bcftools.bcftoolsindex import BcftoolsIndex
        indexer = BcftoolsIndex(camel)
        SnakemakeUtils.add_pickle_inputs(indexer, input)
        step = Step(str(rule), indexer, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(indexer, output)

rule clair3_variant_calling_unzip_vcf:
    """
    Unzips the VCF file.
    """
    input:
        VCF_GZ = rules.clair3_variant_calling_index_vcf_gz.output.VCF_GZ
    output:
        VCF = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_UNFILTERED_VCF
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'unzip_vcf',
        sample_name = config['sample_name']
    run:
        from camel.app.components.filesystemhelper import FileSystemHelper
        from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
        bcftools_view = BcftoolsView(camel)
        SnakemakeUtils.add_pickle_inputs(bcftools_view, input)
        step = Step(str(rule), bcftools_view, camel, params.running_dir, config)
        output_filename = f'variants-{FileSystemHelper.make_valid(params.sample_name)}.vcf'
        bcftools_view.update_parameters(output_format='VCF', compress_output=False, output_filename=output_filename)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_view, output)
