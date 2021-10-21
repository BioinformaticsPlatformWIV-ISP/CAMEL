from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import gene_detection
from camel.scripts.shigellapipeline.snakefile import flexneritype, subspecies_identification


camel = Camel.get_instance()


rule flexneri_call_variants_gtr_promotor:
    """
    Performs variant calling in the gtr promotor region.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        VCF = Path(config['working_dir']) / 'flexneri_type' / 'gtr_promotor' / 'vcf.io',
        BAM = Path(config['working_dir']) / 'flexneri_type' / 'gtr_promotor' / 'bam.io'
    params:
        working_dir = Path(config['working_dir']) / 'flexneri_type' / 'gtr_promotor',
        promotor_fasta = config['flexneri_type']['fasta_gtr_promotor'],
        read_type = config.get('read_type', 'illumina')
    threads: 2
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.components.workflows.variantcallingwrapper import VariantCallingWrapper
        workflow_input = FastqInput.from_fq_dict(Path(input.IO), params.read_type)
        ref_info = {'path': params.promotor_fasta, 'name': 'gtr_promotor'}
        wrapper = VariantCallingWrapper(params.working_dir)
        wrapper.run_workflow(ref_info, 'sample', workflow_input, {'ploidy': 1}, threads)
        SnakemakeUtils.dump_object([wrapper.output.vcf_unfiltered], Path(output.VCF))
        SnakemakeUtils.dump_object([wrapper.output.bam_file], Path(output.BAM))

rule flexneri_call_gtr_promotor_depth:
    """
    Determines the depth of the gtr promotor region.
    """
    input:
        BAM = rules.flexneri_call_variants_gtr_promotor.output.BAM
    output:
        INFORMS = Path(config['working_dir']) / 'flexneri_type' / 'gtr_promotor' / 'informs-depth.io'
    params:
        running_dir = Path(config['working_dir']) / 'flexneri_type' / 'gtr_promotor'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        step = Step(rule, samtools_depth, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)

checkpoint flexneri_type_prepare_reference_files:
    """
    Creates the folder containing the reference files (FASTA, GFF) for the detected loci.
    """
    input:
        VAL_hits = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_ALL_HITS).format(db='flexneri_type')
    output:
        DIR_FASTA = directory(Path(config['working_dir']) / 'flexneri_type' / 'loci')
    params:
        FASTA_ROOT = Path(config['flexneri_type']['fasta_separate'])
    run:
        import logging
        output_dir = Path(output.DIR_FASTA)
        if not output_dir.exists():
            output_dir.mkdir()
        loci = [io.value.locus for io in SnakemakeUtils.load_object(Path(input.VAL_hits))]
        logging.info(f"Hits found for flexneri loci: {loci}")
        fasta_files = []
        dir_by_locus_name = {dir_locus.name: dir_locus for dir_locus in params.FASTA_ROOT.iterdir()}
        for locus in loci:
            output_dir_locus = output_dir / locus
            if not output_dir_locus.exists():
                output_dir_locus.mkdir()
            fasta_ref = dir_by_locus_name[locus] / f'{locus}.fasta'
            SnakemakeUtils.dump_object([ToolIOFile(fasta_ref)], output_dir_locus / 'fasta.io')
            gff = dir_by_locus_name[locus] / f'{locus}.gff'
            SnakemakeUtils.dump_object([ToolIOFile(gff)], output_dir_locus / 'gff.io')

rule flexneri_map_reads:
    """
    Maps the trimmed reads against the reference sequence for the locus. 
    """
    input:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'flexneri_type' / 'loci' / '{flexneri_locus}' / 'fasta.io',
        FASTQ = Path(config['working_dir']) / 'fq_dict.io'
    output:
        SAM = Path(config['working_dir']) / 'flexneri_type' / 'loci' / '{flexneri_locus}' / 'sam.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'flexneri_type' / 'loci' / wildcards.flexneri_locus
    threads: 2
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        bowtie2_map = Bowtie2Map(camel)
        step = Step(rule, bowtie2_map, camel, params.running_dir, config)
        bowtie2_map.update_parameters(threads=threads, no_unal=None, very_sensitive_local=True, sensitive=False, end_to_end=False)
        fasta_as_io_value = [ToolIOValue(io.path) for io in SnakemakeUtils.load_object(Path(input.INDEX_GENOME_PREFIX))]
        bowtie2_map.add_input_files({'INDEX_GENOME_PREFIX': fasta_as_io_value})
        bowtie2_map.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.FASTQ)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bowtie2_map, output)

rule flexneri_sam_to_indexed_bam:
    """
    Converts the read mapping SAM file to an indexed BAM file.
    """
    input:
        SAM = rules.flexneri_map_reads.output.SAM
    output:
        BAM = Path(config['working_dir']) / 'flexneri_type' / 'loci' / '{flexneri_locus}' / 'bam.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'flexneri_type' / 'loci' / wildcards.flexneri_locus
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

rule flexneri_pileup:
    """
    Creates a pileup based on the input BAM file and reference sequence.
    """
    input:
        FASTA = Path(config['working_dir']) / 'flexneri_type' / 'loci' / '{flexneri_locus}' / 'fasta.io',
        BAM = rules.flexneri_sam_to_indexed_bam.output.BAM
    output:
        VCF_GZ = Path(config['working_dir']) / 'flexneri_type' / 'loci' / '{flexneri_locus}' / 'vcf-pileup.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'flexneri_type' / 'loci' / wildcards.flexneri_locus
    run:
        from camel.app.tools.samtools.samtoolsmpileup import SamtoolsMPileup
        pileup = SamtoolsMPileup(camel)
        SnakemakeUtils.add_pickle_inputs(pileup, input)
        step = Step(rule, pileup, camel, params.running_dir, config)
        pileup.update_parameters(output_format='vcf', count_orphans=True)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(pileup, output)

rule flexneri_snp_calling:
    """
    Performs SNP calling for the flexneri loci.
    # TODO: Add Z-score filter?
    """
    input:
        VCF_GZ = rules.flexneri_pileup.output.VCF_GZ
    output:
        VCF_GZ = Path(config['working_dir']) / 'flexneri_type' / 'loci' / '{flexneri_locus}' / 'vcf-gz.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'flexneri_type' / 'loci' / wildcards.flexneri_locus
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

rule flexneri_filter_snps:
    """
    Filters SNPs for the flexneri type.
    """
    input:
        VCF_GZ = rules.flexneri_snp_calling.output.VCF_GZ
    output:
        VCF_GZ = Path(config['working_dir']) / 'flexneri_type' / 'loci' / '{flexneri_locus}' / 'vcf-gz_filtered.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'flexneri_type' / 'loci' / wildcards.flexneri_locus
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

rule flexneri_bcftools_csq:
    """
    Determines the consequence of the detected variants using bcftools csq.
    """
    input:
        VCF_GZ = rules.flexneri_filter_snps.output.VCF_GZ,
        FASTA = Path(config['working_dir']) / 'flexneri_type' / 'loci' / '{flexneri_locus}' / 'fasta.io',
        GFF = Path(config['working_dir']) / 'flexneri_type' / 'loci' / '{flexneri_locus}' / 'gff.io'
    output:
        VCF = Path(config['working_dir']) / 'flexneri_type' / 'loci' / '{flexneri_locus}' / 'vcf-csq.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'flexneri_type' / 'loci' / wildcards.flexneri_locus
    run:
        from camel.app.tools.bcftools.bcftoolscsq import BcftoolsCsq
        bcftools_csq = BcftoolsCsq(camel)
        bcftools_csq.update_parameters(local_csq=None)
        SnakemakeUtils.add_pickle_inputs(bcftools_csq, input)
        step = Step(rule, bcftools_csq, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_csq, output)

rule flexneri_parse_csq:
    """
    Parses the VCF file generated by bcftools csq.
    """
    input:
        VCF = rules.flexneri_bcftools_csq.output.VCF
    output:
        VAL_mut = Path(config['working_dir']) / 'flexneri_type' / 'loci' / '{flexneri_locus}' / 'val-mut.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'flexneri_type' / 'loci' / wildcards.flexneri_locus
    run:
        from camel.app.tools.bcftoolscsqparser.bcftoolscsqparser import CsqParser
        parser = CsqParser(camel)
        SnakemakeUtils.add_pickle_inputs(parser, input)
        step = Step(rule, parser, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(parser, output)

rule flexneri_combine_hits:
    """
    Combines the hits for the different loci.
    """
    input:
        VAL_mut = lambda wildcards: flexneritype.aggregate_input(wildcards, checkpoints, config),
        VCF_csq = lambda wildcards: flexneritype.aggregate_input_vcf(wildcards, checkpoints, config),
        VCF = rules.flexneri_call_variants_gtr_promotor.output.VCF
    output:
        INFORMS = Path(config['working_dir']) / 'flexneri_type' / 'detection' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'flexneri_type',
        dir_fasta = config['flexneri_type']['fasta_combined'],
        tsv_profiles = config['flexneri_type']['tsv_profiles']
    run:
        from camel.app.io.tooliodirectory import ToolIODirectory
        mutations_by_locus = {}
        vcf_by_locus = {}
        for io_path_mut, io_path_csq in zip([Path(x) for x in input.VAL_mut], [Path(x) for x in input.VCF_csq]):
            locus = io_path_mut.parent.name
            mutations_by_locus[locus] = SnakemakeUtils.load_object(io_path_mut)
            vcf_by_locus[locus] = SnakemakeUtils.load_object(io_path_csq)

        from camel.app.tools.pipelines.shigella.flexneritypedetector import FlexneriTypeDetector
        detector = FlexneriTypeDetector(camel)
        detector.add_input_files({f'VAL_mut_{l}': ms for l, ms in mutations_by_locus.items()})
        detector.add_input_files({f'VCF_csq_{l}': ms for l, ms in vcf_by_locus.items()})
        detector.add_input_files({
            'DIR_FASTA': [ToolIODirectory(Path(params.dir_fasta))],
            'TSV': [ToolIOFile(Path(params.tsv_profiles))],
            'VCF': SnakemakeUtils.load_object(Path(input.VCF))
        })
        step = Step(rule, detector, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(detector, output)

rule flexneri_create_report:
    """
    Creates the report section for the flexneri type detection. 
    """
    input:
        INFORMS_subspecies = Path(config['working_dir']) / subspecies_identification.OUTPUT_SPECIES_SUBSPECIES_INFORMS,
        INFORMS_detection = rules.flexneri_combine_hits.output.INFORMS,
        INFORMS_gtr_depth = rules.flexneri_call_gtr_promotor_depth.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / flexneritype.OUTPUT_FLEXNERI_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'flexneri_type' / 'report'
    run:
        from camel.app.tools.pipelines.shigella.flexneritypereporter import FlexneriTypeReporter
        reporter = FlexneriTypeReporter(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule flexneri_report_empty:
    """
    Creates an empty report when (sub)species identification is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / flexneritype.OUTPUT_FLEXNERI_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Flexneri type determination', output.VAL_HTML)

rule flexneri_type_dump_summary_info:
    """
    Dumps the summary information for the flexneri type detection.
    """
    input:
        INFORMS_detection = rules.flexneri_combine_hits.output.INFORMS,
        INFORMS_gtr_depth = rules.flexneri_call_gtr_promotor_depth.output.INFORMS
    output:
        TSV = Path(config['working_dir'], flexneritype.OUTPUT_FLEXNERI_SUMMARY)
    run:
        informs = SnakemakeUtils.load_object(Path(input.INFORMS_detection))
        table_data = [
            ['flexneri_detected_type', informs['detected_type']],
            ['flexneri_loci', informs['loci']]
        ]
        with open(output.TSV, 'w') as handle:
            for key, value in table_data:
                handle.write('{}\t{}'.format(key, str(value)))
                handle.write('\n')
