"""
This Snakefile performs variant calling using the clair3 pipeline.
"""
from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import variant_calling_clair3

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
        SnakemakeUtils.dump_object([ToolIOValue(params.reference['path'])], Path(output.INDEX_GENOME_PREFIX))
        SnakemakeUtils.dump_object([ToolIOFile(Path(params.reference['path']))], Path(output.FASTA))
        SnakemakeUtils.dump_object(params.reference, Path(output.INFORMS))


rule variant_calling_read_mapping:
    """
    Maps the trimmed reads to the assembly.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
        FASTA = rules.variant_calling_prep_reference.output.FASTA
    output:
        SAM = Path(config['working_dir']) / 'variant_calling' / 'read_mapping' / 'sam.io',
        INFORMS = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_MAPPING_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'read_mapping',
        read_type = config.get('read_type', 'illumina')
    threads: 4
    priority: 1
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        minimap2 = Minimap2Mapping(camel)
        step = Step(rule, minimap2, camel, params.running_dir, config)
        minimap2.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(minimap2, 'FASTA', Path(input.FASTA))
        minimap2.add_input_files(SnakePipelineUtils.extracts_fq_input(
            Path(input.IO), key_se='FASTQ', drop_empty=True, read_type='SE'))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(minimap2, output)
        exit()

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
        BAM = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_BAM
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'alignment_sorting'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort(camel)
        step = Step(rule, samtools_sort, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_sort, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_sort, output)

rule variant_calling_index_bam:
    """
    Index the bam file.
    """
    input:
        BAM = rules.variant_calling_alignment_sorting.output.BAM
    output:
        BAM = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_BAM_INDEX
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'alignment_sorting'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        samtools_index = SamtoolsIndex(camel)
        step = Step(rule, samtools_index, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_index, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_index, output)

rule variant_calling_index_fasta:
    """
    Index the bam file.
    """
    input:
        FASTA = rules.variant_calling_prep_reference.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_FASTA_INDEX
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'reference'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools_fastaindex = SamtoolsFastaIndex(camel)
        step = Step(rule, samtools_fastaindex, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_fastaindex, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_fastaindex, output)

rule variant_calling_calculate_depth:
    """
    Calculates the median depth of the alignment.
    """
    input:
        BAM = rules.variant_calling_alignment_sorting.output.BAM
    output:
        TSV = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_DEPTH_TSV,
        INFORMS = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_DEPTH_INFORMS
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

rule clair3_variant_calling:
    """
    This step actually calls variants with clair3
    """
    input:
        FASTA = rules.variant_calling_prep_reference.output.FASTA,
        BAM = rules.variant_calling_alignment_sorting.output.BAM,
        BAM_INDEX = rules.variant_calling_index_bam.output.BAM,
        FASTA_INDEX = rules.variant_calling_index_fasta.output.FASTA
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
        step = Step(rule, clair3, camel, params.running_dir, config)
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

rule variant_calling_normalize_indels:
    """
    Normalizes indels.
    """
    input:
        VCF_GZ = rules.clair3_variant_calling.output.VCF,
        FASTA = rules.variant_calling_prep_reference.output.FASTA
    output:
        VCF_GZ = Path(config['working_dir']) / 'variant_calling' / 'norm' / 'vcf_gz.io',
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'norm'
    run:
        from camel.app.tools.bcftools.bcftoolsnorm import BcftoolsNorm
        bcftools_norm = BcftoolsNorm(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(bcftools_norm, input)
        step = Step(rule, bcftools_norm, camel, params.running_dir, config)
        bcftools_norm.update_parameters(output_format='z')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_norm, output)

rule variant_calling_index_vcf_gz:
    """
    Indexes the VCF file.
    """
    input:
        VCF_GZ = rules.variant_calling_normalize_indels.output.VCF_GZ
    output:
        VCF_GZ = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'norm'
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
        VCF_GZ = rules.variant_calling_index_vcf_gz.output.VCF_GZ
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
        step = Step(rule, bcftools_view, camel, params.running_dir, config)
        output_filename = f'variants-{FileSystemHelper.make_valid(params.sample_name)}.vcf'
        bcftools_view.update_parameters(output_format='VCF', compress_output=False, output_filename=output_filename)
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
        FASTA = Path(config['working_dir']) / variant_calling_clair3.OUTPUT_VARIANT_CALLING_CONSENSUS
    params:
        running_dir = Path(config['working_dir']) / 'variant_calling' / 'consensus',
        output_filename='bcftools_consensus.fasta'
    run:
        from camel.app.tools.bcftools.bcftoolsconsensus import BcftoolsConsensus
        bcftools_consensus = BcftoolsConsensus(camel)
        SnakemakeUtils.add_pickle_inputs(bcftools_consensus, input)
        bcftools_consensus.update_parameters(output_filename = params.output_filename)
        step = Step(rule, bcftools_consensus, camel, Path(params.running_dir), config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_consensus, output)
