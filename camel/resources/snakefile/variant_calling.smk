"""
This Snakefile performs variant calling using the samtools pipeline.
"""

import os

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_READS_PE, OUTPUT_READ_TRIMMING_READS_SE_FWD, \
    OUTPUT_READ_TRIMMING_READS_SE_REV
from camel.resources.snakefile.read_trimming_iontorrent import OUTPUT_TRIMMING_IT_READS
from camel.resources.snakefile.variant_calling import OUTPUT_VARIANT_CALLING_SUMMARY, \
    OUTPUT_VARIANT_CALLING_UNFILTERED_VCF, OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ, OUTPUT_VARIANT_CALLING_REPORT, \
    OUTPUT_VARIANT_CALLING_BAM, OUTPUT_VARIANT_CALLING_CONSENSUS, OUTPUT_VARIANT_CALLING_MAPPING_INFORMS, \
    OUTPUT_VARIANT_CALLING_INFORMS_ALL
from camel.resources.snakefile.variant_filtering import OUTPUT_VARIANT_FILTERING_VCF, OUTPUT_VARIANT_FILTERING_STATS, \
    get_filtering_param

camel = Camel.get_instance()


rule Variant_calling_prep_reference:
    """
    Converts the reference to a Snakemake / CAMEL compatible format. Creates the reference metadata.
    """
    output:
        INDEX_GENOME_PREFIX=os.path.join(config['working_dir'], 'variant_calling', 'reference', 'genome_prefix.io'),
        FASTA=os.path.join(config['working_dir'], 'variant_calling', 'reference', 'fasta.io'),
        INFORMS=os.path.join(config['working_dir'], 'variant_calling', 'reference', 'informs.io')
    params:
        reference=config['variant_calling'].get('reference')
    run:
        SnakemakeUtils.dump_object([ToolIOValue(params.reference['path'])], output.INDEX_GENOME_PREFIX)
        SnakemakeUtils.dump_object([ToolIOFile(params.reference['path'])], output.FASTA)
        SnakemakeUtils.dump_object(params.reference, output.INFORMS)

rule Variant_calling_collect_input:
    """
    Collects the input for the variant calling pipeline.
    """
    input:
        ILLUMINA_FASTQ_PE = os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_PE) if config.get('read_type', 'illumina') == 'illumina' else [],
        ILLUMINA_FASTQ_SE_FWD = os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_SE_FWD) if config.get('read_type', 'illumina') == 'illumina' else [],
        ILLUMINA_FASTQ_SE_REV = os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_SE_REV) if config.get('read_type', 'illumina') == 'illumina' else [],
        IONTORRENT_FASTQ_SE=os.path.join(config['working_dir'], OUTPUT_TRIMMING_IT_READS) if config.get('read_type', 'illumina') == 'iontorrent' else []
    output:
        FASTQ=os.path.join(config['working_dir'], 'variant_calling', 'input-fastq.io')
    params:
        read_type = config.get('read_type', 'illumina')
    run:
        output_dict = {}
        if params.read_type == 'illumina':
            output_dict = {'FASTQ_PE': SnakemakeUtils.load_object(input.ILLUMINA_FASTQ_PE)}
            se_reads = SnakemakeUtils.load_object(input.ILLUMINA_FASTQ_SE_FWD) + \
                       SnakemakeUtils.load_object(input.ILLUMINA_FASTQ_SE_REV)
            if len(se_reads) > 0:
                output_dict['FASTQ_SE'] = se_reads
        else:
            output_dict = {'FASTQ_SE': SnakemakeUtils.load_object(input.IONTORRENT_FASTQ_SE)}
        SnakemakeUtils.dump_object(output_dict, output[0])


rule Variant_calling_read_mapping:
    """
    Maps the trimmed reads to the assembly.
    """
    input:
        FASTQ=os.path.join(config['working_dir'], 'variant_calling', 'input-fastq.io'),
        INDEX_GENOME_PREFIX=os.path.join(config['working_dir'], 'variant_calling', 'reference', 'genome_prefix.io')
    output:
        SAM=os.path.join(config['working_dir'], 'variant_calling', 'read_mapping', 'sam.io'),
        INFORMS=os.path.join(config['working_dir'], OUTPUT_VARIANT_CALLING_MAPPING_INFORMS)
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_calling', 'read_mapping')
    threads: 4
    priority: 1
    run:
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        bowtie2_map = Bowtie2Map(camel)
        step = Step(rule, bowtie2_map, camel, params.running_dir, config)
        bowtie2_map.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(bowtie2_map, 'INDEX_GENOME_PREFIX', input.INDEX_GENOME_PREFIX)
        bowtie2_map.add_input_files(SnakemakeUtils.load_object(input.FASTQ))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bowtie2_map, output)

rule Variant_calling_sam_to_bam:
    """
    Converts the mapped reads SAM file to BAM format.
    """
    input:
        SAM=os.path.join(config['working_dir'], 'variant_calling', 'read_mapping', 'sam.io')
    output:
        BAM=os.path.join(config['working_dir'], 'variant_calling', 'read_mapping', 'bam.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_calling', 'read_mapping')
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        samtools_view = SamtoolsView(camel)
        step = Step(rule, samtools_view, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_view, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_view, output)

rule Variant_calling_alignment_sorting:
    """
    Sorts the alignment.
    """
    input:
        BAM=os.path.join(config['working_dir'], 'variant_calling', 'read_mapping', 'bam.io')
    output:
        BAM=os.path.join(config['working_dir'], OUTPUT_VARIANT_CALLING_BAM)
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_calling', 'alignment_sorting')
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort(camel)
        step = Step(rule, samtools_sort, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_sort, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_sort, output)

rule Variant_calling_mpileup:
    """
    This step creates a multi-way pileup using samtools.
    """
    input:
        FASTA=os.path.join(config['working_dir'], 'variant_calling',  'reference', 'fasta.io'),
        BAM=os.path.join(config['working_dir'], 'variant_calling', 'alignment_sorting', 'bam-sorted.io')
    output:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_calling', 'mpileup', 'vcf_gz.io'),
        INFORMS=os.path.join(config['working_dir'], 'variant_calling', 'mpileup', 'informs.io')
    priority: 1
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_calling', 'mpileup'),
        count_orphans = config['variant_calling'].get('count_orphans', True),
        min_mapping_quality = config['variant_calling'].get('minimal_mq'),
        min_base_quality = config['variant_calling'].get('minimal_bq'),
        disable_baq = config['variant_calling'].get('disable_baq')
    threads: 8
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

rule Variant_calling_bcftools_call:
    """
    Variant calling using bcftools. 
    """
    input:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_calling', 'mpileup', 'vcf_gz.io')
    output:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_calling', 'calling', 'vcf_gz.io'),
        INFORMS=os.path.join(config['working_dir'], 'variant_calling', 'calling', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_calling', 'calling'),
        ploidy=config['variant_calling']['ploidy'],
        calling_method=config['variant_calling'].get('calling_method'),
        skip_variants=config['variant_calling'].get('skip_variants'),
        variants_only=config['variant_calling'].get('variants_only'),
        mutation_rate=config['variant_calling'].get('mutation_rate')
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

rule Variant_calling_VCF_indexing:
    """
    Indexes the VCF file.
    """
    input:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_calling', 'calling', 'vcf_gz.io')
    output:
        VCF_GZ=os.path.join(config['working_dir'], OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ)
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_calling', 'calling')
    run:
        from camel.app.tools.bcftools.bcftoolsindex import BcftoolsIndex
        indexer = BcftoolsIndex(camel)
        SnakemakeUtils.add_pickle_inputs(indexer, input)
        step = Step(rule, indexer, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(indexer, output)

rule Variant_calling_unzip_vcf:
    """
    Unzips the VCF file.
    """
    input:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_calling', 'calling', 'vcf_gz.io')
    output:
        VCF=os.path.join(config['working_dir'], OUTPUT_VARIANT_CALLING_UNFILTERED_VCF)
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_calling', 'unzip_vcf')
    run:
        from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
        bcftools_view = BcftoolsView(camel)
        SnakemakeUtils.add_pickle_inputs(bcftools_view, input)
        step = Step(rule, bcftools_view, camel, params.running_dir, config)
        bcftools_view.update_parameters(output_format='VCF', compress_output=False, output_filename='variants.vcf')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_view, output)

rule Variant_calling_create_consensus:
    """
    Creates the consensus sequence by applying the detected variants to the reference genome.
    """
    input:
        VCF_GZ=os.path.join(config['working_dir'], OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ),
        FASTA=os.path.join(config['working_dir'], 'variant_calling', 'reference', 'fasta.io')
    output:
        FASTA=os.path.join(config['working_dir'], OUTPUT_VARIANT_CALLING_CONSENSUS)
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_calling', 'consensus'),
        output_filename='bcftools_consensus.fasta'
    run:
        from camel.app.tools.bcftools.bcftoolsconsensus import BcftoolsConsensus
        bcftools_consensus = BcftoolsConsensus(camel)
        SnakemakeUtils.add_pickle_inputs(bcftools_consensus, input)
        bcftools_consensus.update_parameters(output_filename=params.output_filename)
        step = Step(rule, bcftools_consensus, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_consensus, output)

rule Variant_calling_report:
    """
    Creates a report for the variant calling.
    """
    input:
        VCF=os.path.join(config['working_dir'], OUTPUT_VARIANT_CALLING_UNFILTERED_VCF),
        VCF_filt=os.path.join(config['working_dir'], OUTPUT_VARIANT_FILTERING_VCF),
        VCF_filt_regions=os.path.join(config['working_dir'], 'variant_filtering', 'regions', 'vcf.io') if get_filtering_param(config, 'region', 'bed_file') is not None else [],
        INFORMS_reference=os.path.join(config['working_dir'], 'variant_calling', 'reference', 'informs.io'),
        INFORMS_mapping=os.path.join(config['working_dir'], 'variant_calling', 'read_mapping', 'informs.io'),
        INFORMS_calling=os.path.join(config['working_dir'], 'variant_calling', 'calling', 'informs.io'),
        JSON=os.path.join(config['working_dir'], OUTPUT_VARIANT_FILTERING_STATS)
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_VARIANT_CALLING_REPORT),
        INFORMS=os.path.join(config['working_dir'], 'variant_calling', 'report', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_calling', 'report'),
        regions_bed_file=get_filtering_param(config, 'region', 'bed_file'),
        sample_name=config['sample_name']
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

rule Variant_calling_dump_summary_info:
    """
    Dumps the summary information from the variant calling workflow.
    """
    input:
        INFORMS_mapping=os.path.join(config['working_dir'], 'variant_calling', 'read_mapping', 'informs.io')
    output:
        os.path.join(config['working_dir'], OUTPUT_VARIANT_CALLING_SUMMARY)
    run:
        informs_mapping = SnakemakeUtils.load_object(input.INFORMS_mapping)
        summary_data = [
            ['vc-mapping_rate', informs_mapping['stats_map_rate']]
        ]
        with open(output[0], 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')

rule Variant_calling_collect_command_informs:
    """
    This rule is used to collect the commands that were used.
    """
    input:
        INFORMS_mapping=os.path.join(config['working_dir'], 'variant_calling', 'read_mapping', 'informs.io'),
        INFORMS_mpileup=os.path.join(config['working_dir'], 'variant_calling', 'mpileup', 'informs.io'),
        INFORMS_calling=os.path.join(config['working_dir'], 'variant_calling', 'calling', 'informs.io')
    output:
        INFORMS_ALL=os.path.join(config['working_dir'], OUTPUT_VARIANT_CALLING_INFORMS_ALL)
    run:
        all_informs = []
        for io_file in input:
            informs = SnakemakeUtils.load_object(io_file)
            informs['_tag'] = 'Variant calling'
            all_informs.append(informs)
        SnakemakeUtils.dump_object(all_informs, output.INFORMS_ALL)
