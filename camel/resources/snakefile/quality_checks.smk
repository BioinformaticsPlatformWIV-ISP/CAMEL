"""
This Snakefile performs a set of advanced quality checks.
"""
import os

from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.assembly_spades import OUTPUT_ASSEMBLY_FASTA
from camel.resources.snakefile.quality_checks import OUTPUT_QUALITY_CHECKS_REPORT, OUTPUT_QUALITY_CHECKS_SUMMARY
from camel.resources.snakefile.read_trimming_iontorrent import get_trimming_folder
from camel.resources.snakefile.sequence_typing import OUTPUT_TYPING_HITS


rule QC_checks_fastQC_additional_checks:
    """
    Tests additional quality metrics based on the FastQC data file output.
    """
    input:
        TXT=os.path.join(get_trimming_folder(config), 'fastqc-post', 'txt.io'),
        TXT_RAW=os.path.join(get_trimming_folder(config), 'fastqc-pre', 'txt.io')
    output:
        INFORMS=os.path.join(config['working_dir'], 'quality_checks', 'fastqc_checks', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'quality_checks', 'fastqc_checks'),
        gc_content=config['quality_checks']['expected_gc_content']
    run:
        from camel.app.tools.fastqc.fastqcadditionalchecks import FastQCAdditionalChecks
        fastqc_checks = FastQCAdditionalChecks(camel)
        step = Step(rule, fastqc_checks, camel, params.running_dir, config)
        fastqc_checks.update_parameters(gc_content_reference=params.gc_content)
        SnakemakeUtils.add_pickle_inputs(fastqc_checks, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastqc_checks, output)

rule QC_checks_bowtie2_index_assembly:
    """
    Creates a Bowtie2 index for the assembly.
    """
    input:
        FASTA_REF=os.path.join(config['working_dir'], OUTPUT_ASSEMBLY_FASTA)
    output:
        INDEX_GENOME_PREFIX=os.path.join(config['working_dir'], 'quality_checks', 'bowtie2_index', 'genome_prefix.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'quality_checks', 'bowtie2_index')
    run:
        from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index
        bowtie2_index = Bowtie2Index(camel)
        step = Step(rule, bowtie2_index, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(bowtie2_index, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bowtie2_index, output)

rule QC_checks_map_against_assembly:
    """
    Maps the trimmed reads to the assembly.
    """
    input:
        FASTQ=os.path.join(config['working_dir'], 'variant_calling', 'input-fastq.io'),
        INDEX_GENOME_PREFIX=os.path.join(config['working_dir'], 'quality_checks', 'bowtie2_index', 'genome_prefix.io')
    output:
        SAM=os.path.join(config['working_dir'], 'quality_checks', 'read_mapping', 'sam.io'),
        INFORMS=os.path.join(config['working_dir'], 'quality_checks', 'read_mapping', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'quality_checks', 'read_mapping')
    priority: 1
    run:
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        bowtie2_map = Bowtie2Map(camel)
        step = Step(rule, bowtie2_map, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_input(bowtie2_map, 'INDEX_GENOME_PREFIX', input.INDEX_GENOME_PREFIX)
        bowtie2_map.add_input_files(SnakemakeUtils.load_object(input.FASTQ))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bowtie2_map, output)

rule QC_checks_SAM_to_BAM:
    """
    Converts the mapped reads SAM file to BAM format.
    """
    input:
        SAM=os.path.join(config['working_dir'], 'quality_checks', 'read_mapping', 'sam.io')
    output:
        BAM=os.path.join(config['working_dir'], 'quality_checks', 'read_mapping', 'bam.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'quality_checks', 'read_mapping')
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        samtools_view = SamtoolsView(camel)
        step = Step(rule, samtools_view, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_view, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_view, output)

rule QC_checks_sort_BAM:
    """
    Sorts the alignment.
    """
    input:
        BAM=os.path.join(config['working_dir'], 'quality_checks', 'read_mapping', 'bam.io')
    output:
        BAM=os.path.join(config['working_dir'], 'quality_checks', 'read_mapping', 'bam-sorted.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'quality_checks', 'read_mapping')
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort(camel)
        step = Step(rule, samtools_sort, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_sort, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_sort, output)

rule QC_checks_calculate_coverage:
    """
    Determines the read depth at every position to estimate the coverage.
    """
    input:
        BAM=os.path.join(config['working_dir'], 'quality_checks', 'read_mapping', 'bam-sorted.io')
    output:
        INFORMS=os.path.join(config['working_dir'], 'quality_checks', 'coverage_calculation', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'quality_checks', 'coverage_calculation')
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth(camel)
        step = Step(rule, samtools_depth, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)

rule QC_checks_get_cgMLST_stats:
    """
    Retrieves the number of cgMLST genes that were detected (only perfect hits are considered).
    """
    input:
        hits=os.path.join(config['working_dir'], OUTPUT_TYPING_HITS.format(
            locus_type='DNA',
            scheme=config['quality_checks']['typing_scheme'],
            detection_method=config['detection_method'])),
        INFORMS_scheme=os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-locus_set.io').format(
            scheme=config['quality_checks']['typing_scheme'])
    output:
        INFORMS=os.path.join(config['working_dir'], 'quality_checks', 'metrics', 'informs-st.io')
    params:
        scheme_name=config['quality_checks']['typing_scheme']
    run:
        title = SnakemakeUtils.load_object(input.INFORMS_scheme)['title']
        all_hits = SnakemakeUtils.load_object(input.hits)
        nb_perfect = len([v for v in all_hits if v.value.is_perfect_hit()])
        SnakemakeUtils.dump_object(
            {'hits_found': nb_perfect,
             'nb_of_loci': len(all_hits),
             'title': title,
             'scheme_name': params.scheme_name},
            output.INFORMS)

rule QC_checks_check_additional_metrics:
    """
    Checks if the quality criteria were met.
    """
    input:
        INFORMS_coverage=os.path.join(config['working_dir'], 'quality_checks', 'coverage_calculation', 'informs.io'),
        INFORMS_fastqc_checks=os.path.join(config['working_dir'], 'quality_checks', 'fastqc_checks', 'informs.io'),
        INFORMS_mapping=os.path.join(config['working_dir'], 'quality_checks', 'read_mapping', 'informs.io'),
        INFORMS_cgmlst=os.path.join(config['working_dir'], 'quality_checks', 'metrics', 'informs-st.io')
    output:
        INFORMS=os.path.join(config['working_dir'], 'quality_checks', 'metrics', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'quality_checks')
    run:
        from camel.app.tools.pipelines.quality_checks.qualitycriteriachecker import QualityCriteriaChecker
        quality_checker = QualityCriteriaChecker(camel)
        step = Step(rule, quality_checker, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(quality_checker, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quality_checker, output)

rule QC_checks_visualize_qc_checks:
    """
    Creates plots to visualize the QC checks, shows the warning and fail threshold along with the detected value.
    """
    input:
        INFORMS_cov=os.path.join(config['working_dir'], 'quality_checks', 'coverage_calculation', 'informs.io'),
        INFORMS_mlst = os.path.join(config['working_dir'], 'quality_checks', 'metrics', 'informs-st.io'),
        INFORMS_map = os.path.join(config['working_dir'], 'quality_checks', 'read_mapping', 'informs.io'),
    output:
        PNG_cov=os.path.join(config['working_dir'], 'quality_checks', 'visualization', 'png-cov.io'),
        PNG_st=os.path.join(config['working_dir'], 'quality_checks', 'visualization', 'png-st.io'),
        PNG_mapping=os.path.join(config['working_dir'], 'quality_checks', 'visualization', 'png-mapping.io'),
    params:
        running_dir=os.path.join(config['working_dir'], 'quality_checks', 'visualization')
    run:
        from camel.app.tools.pipelines.quality_checks.qcvisualization import QCVisualization
        visuzalization = QCVisualization(camel)
        SnakemakeUtils.add_pickle_inputs(visuzalization, input)
        step = Step(rule, visuzalization, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(visuzalization, output)

rule QC_checks_report:
    """
    Creates a report containing the quality checks information.
    """
    input:
        PNG_cov=os.path.join(config['working_dir'], 'quality_checks', 'visualization', 'png-cov.io'),
        PNG_st=os.path.join(config['working_dir'], 'quality_checks', 'visualization', 'png-st.io'),
        PNG_mapping=os.path.join(config['working_dir'], 'quality_checks', 'visualization', 'png-mapping.io'),
        INFORMS_coverage=os.path.join(config['working_dir'], 'quality_checks', 'coverage_calculation', 'informs.io'),
        INFORMS_mapping=os.path.join(config['working_dir'], 'quality_checks', 'read_mapping', 'informs.io'),
        INFORMS_additional_checks=os.path.join(config['working_dir'], 'quality_checks', 'fastqc_checks', 'informs.io'),
        INFORMS_quality_criteria=os.path.join(config['working_dir'], 'quality_checks', 'metrics', 'informs.io'),
        INFORMS_st_genes=os.path.join(config['working_dir'], 'quality_checks', 'metrics', 'informs-st.io'),
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_QUALITY_CHECKS_REPORT)
    params:
        running_dir=os.path.join(config['working_dir'], 'quality_checks', 'report'),
    run:
        from camel.app.tools.pipelines.quality_checks.htmlreporterqualitychecks import HtmlReporterQualityChecks
        reporter = HtmlReporterQualityChecks(camel)
        step = Step(rule, reporter, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule QC_checks_dump_summary_info:
    """ 
    Dumps the summary information in tabular format.
    """
    input:
        INFORMS_coverage=os.path.join(config['working_dir'], 'quality_checks', 'coverage_calculation', 'informs.io'),
        INFORMS_mapping=os.path.join(config['working_dir'], 'quality_checks', 'read_mapping', 'informs.io'),
        INFORMS_cgmlst=os.path.join(config['working_dir'], 'quality_checks', 'metrics', 'informs-st.io'),
        INFORMS_additional_checks=os.path.join(config['working_dir'], 'quality_checks', 'fastqc_checks', 'informs.io'),
    output:
        os.path.join(config['working_dir'], OUTPUT_QUALITY_CHECKS_SUMMARY)
    params:
        read_type=config.get('read_type', 'illumina')
    run:
        informs_cgmlst = SnakemakeUtils.load_object(input.INFORMS_cgmlst)
        informs_qc_checks = SnakemakeUtils.load_object(input.INFORMS_additional_checks)['tests']
        summary_data = [
            ('coverage', SnakemakeUtils.load_object(input.INFORMS_coverage)['median_depth']),
            ('assembly_mapping_rate', SnakemakeUtils.load_object(input.INFORMS_mapping)['stats_map_rate']),
            ('st_hits_found', informs_cgmlst['hits_found']),
            ('st_nb_loci', informs_cgmlst['nb_of_loci'])
        ]

        # Add FastQC info
        check_keys = [('check_avg_qs', 'Average quality score'), ('check_gc', 'GC content'),
                      ('check_n_fraction', 'Maximal N-fraction'), ('check_qs_drop', 'Mean Q-score drop'),
                      ('check_slen_dist', 'Sequence length distribution')]
        read_keys = ['_fwd', '_rev'] if params.read_type == 'illumina' else ['']
        for summary_key, inform_key in check_keys:
            for i in range(0, len(read_keys)):
                summary_data.append(('{}{}'.format(summary_key, read_keys[i]), informs_qc_checks[inform_key][i]))
        with open(output[0], 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')
