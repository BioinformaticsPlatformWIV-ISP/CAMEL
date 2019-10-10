from typing import List, Dict

import os

from camel.app.io.tooliovalue import ToolIOValue
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.gene_detection import OUTPUT_GENE_DETECTION_ALL_HITS
from camel.scripts.shigellapipeline.snakefile.flexneritype import OUTPUT_FLEXNERI_SUMMARY, OUTPUT_FLEXNERI_REPORT, \
    OUTPUT_FLEXNERI_REPORT_EMPTY

rule Flexneri_call_variants_gtr_promotor:
    """
    Performs variant calling in the gtr promotor region.
    # TODO, add filtering?
    """
    input:
        FASTQ = os.path.join(config['working_dir'], 'variant_calling', 'input-fastq.io')
    output:
        VCF = os.path.join(config['working_dir'], 'flexneri_type', 'gtr_promotor', 'vcf.io'),
        BAM = os.path.join(config['working_dir'], 'flexneri_type', 'gtr_promotor', 'bam.io')
    params:
        working_dir = os.path.join(config['working_dir'], 'flexneri_type', 'gtr_promotor'),
        promotor_fasta = config['flexneri_type']['fasta_gtr_promotor']
    threads: 2
    run:
        from camel.app.components.workflows.variantcallingwrapper import VariantCallingWrapper
        from camel.app.components.workflows.variantcallingwrapper import VariantCallingInput
        fastq_files = SnakemakeUtils.load_object(input.FASTQ)
        workflow_input = VariantCallingInput(
            pe_reads=fastq_files['FASTQ_PE'],
            se_reads_fwd=fastq_files['FASTQ_SE'][0],
            se_reads_rev=fastq_files['FASTQ_SE'][1]
        )
        ref_info = {'path': params.promotor_fasta, 'name': 'gtr_promotor'}
        wrapper = VariantCallingWrapper(params.working_dir)
        wrapper.run_workflow(ref_info, 'sample', workflow_input, {'ploidy': 1}, threads)
        SnakemakeUtils.dump_object([wrapper.output.vcf_unfiltered], output.VCF)
        SnakemakeUtils.dump_object([wrapper.output.bam_file], output.BAM)

rule Flexneri_call_gtr_promotor_depth:
    """
    Determines the depth of the gtr promotor region.
    """
    input:
        BAM = os.path.join(config['working_dir'], 'flexneri_type', 'gtr_promotor', 'bam.io')
    output:
        INFORMS = os.path.join(config['working_dir'], 'flexneri_type', 'gtr_promotor', 'informs-depth.io')
    params:
        running_dir = os.path.join(config['working_dir'], 'flexneri_type', 'gtr_promotor')
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        step = Step(rule, samtools_depth, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)

checkpoint Flexneri_type_prepare_reference_files:
    """
    Creates the folder containing the reference files (FASTA, GFF) for the detected loci.
    """
    input:
        VAL_hits = os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_ALL_HITS.format(db='flexneri_type'))
    output:
        DIR_FASTA = directory(os.path.join(config['working_dir'], 'flexneri_type', 'loci'))
    params:
        FASTA_ROOT = config['flexneri_type']['fasta_separate']
    run:
        import logging
        from camel.app.io.tooliofile import ToolIOFile
        os.makedirs(output.DIR_FASTA)
        loci = [io.value.locus for io in SnakemakeUtils.load_object(input.VAL_hits)]
        logging.info(f"Hits found for flexneri loci: {loci}")
        fasta_files = []
        dir_by_locus_name = {locus: os.path.join(params.FASTA_ROOT, locus) for locus in os.listdir(params.FASTA_ROOT)}
        for locus in loci:
            output_dir = os.path.join(output.DIR_FASTA, locus)
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir)
            fasta_ref = ToolIOFile(os.path.join(dir_by_locus_name[locus], f'{locus}.fasta'))
            SnakemakeUtils.dump_object([fasta_ref], os.path.join(output_dir, 'fasta.io'))
            gff = ToolIOFile(os.path.join(dir_by_locus_name[locus], f'{locus}.gff'))
            SnakemakeUtils.dump_object([gff], os.path.join(output_dir, 'gff.io'))

rule Flexneri_map_reads:
    """
    Maps the trimmed reads against the reference sequence for the locus. 
    """
    input:
        INDEX_GENOME_PREFIX = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'fasta.io'),
        FASTQ = os.path.join(config['working_dir'], 'variant_calling', 'input-fastq.io')
    output:
        SAM = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'sam.io')
    params:
        running_dir = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}')
    threads: 2
    run:
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        bowtie2_map = Bowtie2Map(camel)
        step = Step(rule, bowtie2_map, camel, params.running_dir, config)
        bowtie2_map.update_parameters(threads=threads, no_unal=None, very_sensitive_local=True, sensitive=False, end_to_end=False)
        fasta_as_io_value = [ToolIOValue(io.path) for io in SnakemakeUtils.load_object(input.INDEX_GENOME_PREFIX)]
        bowtie2_map.add_input_files({'INDEX_GENOME_PREFIX': fasta_as_io_value})
        bowtie2_map.add_input_files(SnakemakeUtils.load_object(input.FASTQ))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bowtie2_map, output)

rule Flexneri_sam_to_indexed_bam:
    """
    Converts the read mapping SAM file to an indexed BAM file.
    """
    input:
        SAM = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'sam.io')
    output:
        BAM = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'bam.io')
    params:
        running_dir = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}')
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        # Convert to BAM
        samtools_view = SamtoolsView(camel)
        step = Step(rule, samtools_view, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_view, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_view, output)

        # Sort BAM
        samtools_sort = SamtoolsSort(camel)
        step = Step(rule, samtools_sort, camel, params.running_dir, config)
        samtools_sort.add_input_files({'BAM': samtools_view.tool_outputs['BAM']})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_sort, output)

rule Flexneri_pileup:
    """
    Creates a pileup based on the input BAM file and reference sequence.
    """
    input:
        FASTA = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'fasta.io'),
        BAM = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'bam.io')
    output:
        VCF_GZ = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'vcf-pileup.io')
    params:
        running_dir = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}')
    run:
        from camel.app.tools.samtools.samtoolsmpileup import SamtoolsMPileup
        pileup = SamtoolsMPileup(camel)
        SnakemakeUtils.add_pickle_inputs(pileup, input)
        step = Step(rule, pileup, camel, params.running_dir, config)
        pileup.update_parameters(output_format='vcf', count_orphans=True)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(pileup, output)

rule Flexneri_snp_calling:
    """
    Performs SNP calling for the flexneri loci.
    # TODO: Add Z-score filter?
    """
    input:
        VCF_GZ = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'vcf-pileup.io')
    output:
        VCF_GZ = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'vcf-gz.io')
    params:
        running_dir = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}')
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
            ploidy=1
        )
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(variant_caller, output)

rule Flexneri_filter_snps:
    input:
        VCF_GZ = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'vcf-gz.io')
    output:
        VCF_GZ = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'vcf-gz_filtered.io')
    params:
        running_dir = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}')
    run:
        from camel.app.tools.variantfiltering.depthfilter import DepthFilter
        depth_filter =  DepthFilter(camel)
        SnakemakeUtils.add_pickle_inputs(depth_filter, input)
        depth_filter.update_parameters(min_depth=10, min_forward_depth=1, min_reverse_depth=1)
        step = Step(rule, depth_filter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(depth_filter, output)

        from camel.app.tools.variantfiltering.mappingqualityfilter import MappingQualityFilter
        mapping_filter = MappingQualityFilter(camel)
        mapping_filter.add_input_files({'VCF_GZ': depth_filter.tool_outputs['VCF_GZ']})
        mapping_filter.run(params.running_dir)

        from camel.app.tools.variantfiltering.snpqualityfilter import SnpQualityFilter
        qual_filter = SnpQualityFilter(camel)
        qual_filter.add_input_files({'VCF_GZ': mapping_filter.tool_outputs['VCF_GZ']})
        qual_filter.run(params.running_dir)
        SnakemakeUtils.dump_tool_outputs(qual_filter, output)

rule Flexneri_bcftools_csq:
    """
    Determines the consequence of the detected variants using bcftools csq.
    """
    input:
        VCF_GZ = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'vcf-gz_filtered.io'),
        # VCF_GZ = os.path.join(config['working_dir'], 'flexneri_type', '{flexneri_locus}', 'vcf-gz.io'),
        FASTA = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'fasta.io'),
        GFF = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'gff.io')
    output:
        VCF = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'vcf-csq.io')
    params:
        running_dir = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}')
    run:
        from camel.app.tools.bcftools.bcftoolscsq import BcftoolsCsq
        bcftools_csq = BcftoolsCsq(camel)
        bcftools_csq.update_parameters(local_csq=None)
        SnakemakeUtils.add_pickle_inputs(bcftools_csq, input)
        step = Step(rule, bcftools_csq, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_csq, output)

rule Flexneri_parse_csq:
    """
    Parses the VCF file generated by bcftools csq.
    """
    input:
        VCF = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'vcf-csq.io')
    output:
        VAL_mut = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'val-mut.io')
    params:
        running_dir = os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}')
    run:
        from camel.app.tools.bcftoolscsqparser.bcftoolscsqparser import CsqParser
        parser = CsqParser(camel)
        SnakemakeUtils.add_pickle_inputs(parser, input)
        step = Step(rule, parser, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(parser, output)

def aggregate_input(wildcards: Dict) -> List[str]:
    """
    Aggregates the input from the read mapping based on the detected loci.
    :param wildcards: Wildcards
    :return: List of input files
    """
    dir_fasta = checkpoints.Flexneri_type_prepare_reference_files.get(**wildcards).output.DIR_FASTA
    loci = [l for l in os.listdir(dir_fasta) if not l.startswith('.')]
    return expand(os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'val-mut.io'), flexneri_locus = loci)

def aggregate_input_vcf(wildcards: Dict) -> List[str]:
    """
    Aggregates the input from the read mapping based on the detected loci.
    :param wildcards: Wildcards
    :return: List of input files
    """
    dir_fasta = checkpoints.Flexneri_type_prepare_reference_files.get(**wildcards).output.DIR_FASTA
    loci = [l for l in os.listdir(dir_fasta) if not l.startswith('.')]
    return expand(os.path.join(config['working_dir'], 'flexneri_type', 'loci', '{flexneri_locus}', 'vcf-csq.io'),
                  flexneri_locus=loci)

rule Flexneri_combine_hits:
    """
    
    """
    input:
        VAL_mut = aggregate_input,
        VCF_csq = aggregate_input_vcf,
        VCF = os.path.join(config['working_dir'], 'flexneri_type', 'gtr_promotor', 'vcf.io')
    output:
        INFORMS = os.path.join(config['working_dir'], 'flexneri_type', 'detection', 'informs.io')
    params:
        running_dir = os.path.join(config['working_dir'], 'flexneri_type'),
        dir_fasta = config['flexneri_type']['fasta_combined'],
        tsv_profiles = config['flexneri_type']['tsv_profiles']
    run:
        from camel.app.io.tooliodirectory import ToolIODirectory
        mutations_by_locus = {}
        vcf_by_locus = {}
        for io_path_mut, io_path_csq in zip(input.VAL_mut, input.VCF_csq):
            locus = os.path.basename(os.path.dirname(io_path_mut))
            mutations_by_locus[locus] = SnakemakeUtils.load_object(io_path_mut)
            vcf_by_locus[locus] = SnakemakeUtils.load_object(io_path_csq)

        from camel.app.tools.pipelines.shigella.flexneritypedetector import FlexneriTypeDetector
        detector = FlexneriTypeDetector(camel)
        detector.add_input_files({f'VAL_mut_{l}': ms for l, ms in mutations_by_locus.items()})
        detector.add_input_files({f'VCF_csq_{l}': ms for l, ms in vcf_by_locus.items()})
        detector.add_input_files({
            'DIR_FASTA': [ToolIODirectory(params.dir_fasta)],
            'TSV': [ToolIOFile(params.tsv_profiles)],
            'VCF': SnakemakeUtils.load_object(input.VCF)
        })
        step = Step(rule, detector, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(detector, output)

        # from camel.app.tools.pipelines.shigella.flexneritypereporter import FlexneriTypeReporter
        # reporter = FlexneriTypeReporter(camel)
        # reporter.add_input_files({f'VAL_mut_{l}': ms for l, ms in mutations_by_locus.items()})
        # reporter.run(params.running_dir)

rule Flexneri_create_report:
    """
    Creates the report section for the flexneri type detection. 
    """
    input:
        INFORMS_subspecies = os.path.join(config['working_dir'], 'subspecies_identification', 'subspecies', 'informs.io'),
        INFORMS_detection = os.path.join(config['working_dir'], 'flexneri_type', 'detection', 'informs.io'),
        INFORMS_gtr_depth = os.path.join(config['working_dir'], 'flexneri_type', 'gtr_promotor', 'informs-depth.io')
    output:
        VAL_HTML = os.path.join(config['working_dir'], OUTPUT_FLEXNERI_REPORT)
    params:
        running_dir = os.path.join(config['working_dir'], 'flexneri_type', 'report')
    run:
        from camel.app.tools.pipelines.shigella.flexneritypereporter import FlexneriTypeReporter
        reporter = FlexneriTypeReporter(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule Flexneri_report_empty:
    """
    Creates an empty report when (sub)species identification is disabled.
    """
    output:
        VAL_HTML = os.path.join(config['working_dir'], OUTPUT_FLEXNERI_REPORT_EMPTY)
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Flexneri type determination', output.VAL_HTML)

rule Flexneri_type_dump_summary_info:
    """
    Dumps the summary information for the flexneri type detection.
    """
    input:
        INFORMS_detection = os.path.join(config['working_dir'], 'flexneri_type', 'detection', 'informs.io'),
        INFORMS_gtr_depth = os.path.join(config['working_dir'], 'flexneri_type', 'gtr_promotor', 'informs-depth.io')
    output:
        os.path.join(config['working_dir'], OUTPUT_FLEXNERI_SUMMARY)
    run:
        informs = SnakemakeUtils.load_object(input.INFORMS_detection)
        table_data = [
            ['flexneri_detected_type', informs['detected_type']],
            ['flexneri_loci', informs['loci']]
        ]
        with open(output[0], 'w') as handle:
            for key, value in table_data:
                handle.write('{}\t{}'.format(key, str(value)))
                handle.write('\n')
