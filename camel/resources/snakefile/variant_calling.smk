"""
This Snakefile performs variant calling using the samtools pipeline.
"""
from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import variant_calling, variant_filtering

camel = Camel.get_instance()


rule variant_calling_prep_reference:
    """
    Converts the reference to a Snakemake / CAMEL compatible format. Creates the reference metadata.
    """
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'variant_calling' / 'reference' / 'genome_prefix.io',
        FASTA = Path(config['working_dir']) / 'variant_calling' / 'reference' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / 'variant_calling' / 'reference' / 'informs.io'
    params:
        reference = config['variant_calling'].get('reference')
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.io.tooliofile import ToolIOFile
        SnakemakeUtils.dump_object([ToolIOValue(params.reference['path'])], output.INDEX_GENOME_PREFIX)
        SnakemakeUtils.dump_object([ToolIOFile(params.reference['path'])], output.FASTA)
        SnakemakeUtils.dump_object(params.reference, output.INFORMS)


rule variant_calling_read_mapping:
    """
    Maps the trimmed reads to the assembly.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
        INDEX_GENOME_PREFIX = rules.variant_calling_prep_reference.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'variant_calling' / 'read_mapping' / 'sam.io',
        INFORMS = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_MAPPING_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'read_mapping'
    threads: 4
    priority: 1
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        bowtie2_map = Bowtie2Map(camel)
        step = Step(rule, bowtie2_map, camel, params.running_dir, config)
        bowtie2_map.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(bowtie2_map, 'INDEX_GENOME_PREFIX', input.INDEX_GENOME_PREFIX)
        bowtie2_map.add_input_files(SnakePipelineUtils.extracts_fq_input(input.IO))
        import pprint
        pprint.pprint(bowtie2_map._tool_inputs)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bowtie2_map, output)

rule variant_calling_sam_to_bam:
    """
    Converts the mapped reads SAM file to BAM format.
    """
    input:
        SAM = rules.variant_calling_read_mapping.output.SAM
    output:
        BAM = Path(config['working_dir']) / 'variant_calling' / 'read_mapping' / 'bam.io'
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'read_mapping'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        samtools_view = SamtoolsView(camel)
        step = Step(rule, samtools_view, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_view, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_view, output)

rule variant_calling_alignment_sorting:
    """
    Sorts the alignment.
    """
    input:
        BAM = rules.variant_calling_sam_to_bam.output.BAM
    output:
        BAM = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_BAM
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'alignment_sorting'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort(camel)
        step = Step(rule, samtools_sort, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_sort, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_sort, output)

rule variant_calling_calculate_depth:
    """
    Calculates the median depth of the alignment.
    """
    input:
        BAM = rules.variant_calling_alignment_sorting.output.BAM
    output:
        INFORMS = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_DEPTH_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'depth'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth(camel)
        step = Step(rule, samtools_depth, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        samtools_depth.update_parameters(output_all_positions_absolutely=None)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)

rule variant_calling_mpileup:
    """
    This step creates a multi-way pileup using samtools.
    """
    input:
        FASTA = rules.variant_calling_prep_reference.output.FASTA,
        BAM = rules.variant_calling_alignment_sorting.output.BAM
    output:
        VCF_GZ = Path(config['working_dir']) / 'variant_calling' / 'mpileup' / 'vcf_gz.io',
        INFORMS = Path(config['working_dir']) / 'variant_calling' / 'mpileup' / 'informs.io'
    priority: 1
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'mpileup',
        count_orphans = config['variant_calling'].get('count_orphans', True),
        min_mapping_quality = config['variant_calling'].get('minimal_mq'),
        min_base_quality = config['variant_calling'].get('minimal_bq'),
        disable_baq = config['variant_calling'].get('disable_baq')
    threads: 1
    run:
        from camel.app.tools.samtools.samtoolsmpileup import SamtoolsMPileup
        samtools_mpileup = SamtoolsMPileup(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_mpileup, input)
        step = Step(rule, samtools_mpileup, camel, params.running_dir, config)
        samtools_mpileup.update_parameters(output_format='vcf')
        if params.count_orphans is not None:
            samtools_mpileup.update_parameters(count_orphans=params.count_orphans)
        if params.min_mapping_quality is not None:
            samtools_mpileup.update_parameters(count_orphans=params.min_mapping_quality)
        if params.min_base_quality is not None:
            samtools_mpileup.update_parameters(count_orphans=params.min_base_quality)
        if params.disable_baq is not None:
            samtools_mpileup.update_parameters(count_orphans=params.disable_baq)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_mpileup, output)

rule variant_calling_bcftools_call:
    """
    Variant calling using bcftools. 
    """
    input:
        VCF_GZ = rules.variant_calling_mpileup.output.VCF_GZ
    output:
        VCF_GZ = Path(config['working_dir']) / 'variant_calling' / 'calling' / 'vcf_gz.io',
        INFORMS = Path(config['working_dir']) / 'variant_calling' / 'calling' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'calling',
        ploidy = config['variant_calling'].get('ploidy', 1),
        calling_method = config['variant_calling'].get('calling_method'),
        skip_variants = config['variant_calling'].get('skip_variants'),
        variants_only = config['variant_calling'].get('variants_only'),
        mutation_rate = config['variant_calling'].get('mutation_rate')
    run:
        from camel.app.tools.bcftools.bcftoolscall import BcftoolsCall
        variant_caller = BcftoolsCall(camel)
        SnakemakeUtils.add_pickle_inputs(variant_caller, input)
        step = Step(rule, variant_caller, camel, params.running_dir, config)
        variant_caller.update_parameters(
            output_format='VCF',
            output_filename='variants.vcf.gz',
            variants_only=True,
            compress_output=True,
            ploidy=params.ploidy
        )
        if params.calling_method is not None:
            variant_caller.update_parameters(calling_method=params.calling_method)
        if params.skip_variants is not None:
            variant_caller.update_parameters(skip_variants=params.skip_variants)
        if params.variants_only is not None:
            variant_caller.update_parameters(variants_only=params.variants_only)
        if params.mutation_rate is not None:
            variant_caller.update_parameters(mutation_rate=params.mutation_rate)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(variant_caller, output)

rule variant_calling_index_vcf_gz:
    """
    Indexes the VCF file.
    """
    input:
        VCF_GZ = rules.variant_calling_bcftools_call.output.VCF_GZ
    output:
        VCF_GZ = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'calling'
    run:
        from camel.app.tools.bcftools.bcftoolsindex import BcftoolsIndex
        indexer = BcftoolsIndex(camel)
        SnakemakeUtils.add_pickle_inputs(indexer, input)
        step = Step(rule, indexer, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(indexer, output)

rule variant_calling_unzip_vcf:
    """
    Unzips the VCF file.
    """
    input:
        VCF_GZ = rules.variant_calling_bcftools_call.output.VCF_GZ
    output:
        VCF = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_UNFILTERED_VCF
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'unzip_vcf'
    run:
        from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
        bcftools_view = BcftoolsView(camel)
        SnakemakeUtils.add_pickle_inputs(bcftools_view, input)
        step = Step(rule, bcftools_view, camel, params.running_dir, config)
        bcftools_view.update_parameters(output_format='VCF', compress_output=False, output_filename='variants.vcf')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_view, output)

rule variant_calling_create_consensus:
    """
    Creates the consensus sequence by applying the detected variants to the reference genome.
    """
    input:
        VCF_GZ = rules.variant_calling_index_vcf_gz.output.VCF_GZ,
        FASTA = rules.variant_calling_prep_reference.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_CONSENSUS
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'consensus',
        output_filename='bcftools_consensus.fasta'
    run:
        from camel.app.tools.bcftools.bcftoolsconsensus import BcftoolsConsensus
        bcftools_consensus = BcftoolsConsensus(camel)
        SnakemakeUtils.add_pickle_inputs(bcftools_consensus, input)
        bcftools_consensus.update_parameters(output_filename = params.output_filename)
        step = Step(rule, bcftools_consensus, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_consensus, output)

rule variant_calling_report:
    """
    Creates a report for the variant calling.
    """
    input:
        VCF = rules.variant_calling_unzip_vcf.output.VCF,
        VCF_filt = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_VCF,
        VCF_filt_regions = Path(config['working_dir'], 'variant_filtering', 'regions', 'vcf.io') if variant_filtering.get_filtering_param(config, 'region', 'bed_file') is not None else [],
        INFORMS_reference = rules.variant_calling_prep_reference.output.INFORMS,
        INFORMS_mapping = rules.variant_calling_read_mapping.output.INFORMS,
        INFORMS_calling = rules.variant_calling_bcftools_call.output.INFORMS,
        INFORMS_depth = rules.variant_calling_calculate_depth.output.INFORMS,
        JSON = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_STATS
    output:
        VAL_HTML = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_REPORT,
        INFORMS = Path(config['working_dir']) / 'variant_calling' / 'report' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'report',
        regions_bed_file = variant_filtering.get_filtering_param(config, 'region', 'bed_file'),
        sample_name = config['sample_name']
    run:
        from camel.app.tools.pipelines.variant_calling.variantcallingreporter import VariantCallingReporter
        reporter = VariantCallingReporter(camel)
        step = Step(rule, reporter, camel, params.running_dir, config)
        keys = [k for k in input.keys() if k != 'VCF_filt_regions']
        SnakemakeUtils.add_pickle_inputs(reporter, input, keys=keys)
        if params.regions_bed_file is not None:
            reporter.add_input_files({'BED': [ToolIOFile(params.regions_bed_file)]})
            SnakemakeUtils.add_pickle_input(reporter, 'VCF_filt_regions', input.VCF_filt_regions)
        reporter.add_input_files({'VAL_Sample': [ToolIOValue(params.sample_name)]})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule variant_calling_dump_summary_info:
    """
    Dumps the summary information from the variant calling workflow.
    """
    input:
        INFORMS_mapping = rules.variant_calling_read_mapping.output.INFORMS,
        INFORMS_depth = rules.variant_calling_calculate_depth.output.INFORMS
    output:
        Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_SUMMARY
    run:
        informs_mapping = SnakemakeUtils.load_object(input.INFORMS_mapping)
        informs_depth = SnakemakeUtils.load_object(input.INFORMS_depth)
        summary_data = [
            ['vc-mapping_rate', informs_mapping['stats_map_rate']],
            ['vc-median_depth', informs_depth['median_depth']]
        ]
        with open(output[0], 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')

rule variant_calling_collect_command_informs:
    """
    This rule is used to collect the commands that were used.
    """
    input:
        INFORMS_mapping = rules.variant_calling_read_mapping.output.INFORMS,
        INFORMS_mpileup = rules.variant_calling_mpileup.output.INFORMS,
        INFORMS_calling = rules.variant_calling_bcftools_call.output.INFORMS
    output:
        INFORMS_ALL = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_INFORMS_ALL
    run:
        all_informs = []
        for io_file in input:
            informs = SnakemakeUtils.load_object(io_file)
            informs['_tag'] = 'Variant calling'
            all_informs.append(informs)
        SnakemakeUtils.dump_object(all_informs, output.INFORMS_ALL)
