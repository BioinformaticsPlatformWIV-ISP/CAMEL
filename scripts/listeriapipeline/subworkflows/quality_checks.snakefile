QUALITY_CHECKS_WORKING_DIR=os.path.join(__WORKING_DIR, 'quality_checks')
QUALITY_CHECKS_REPORT=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'report-html.io')
QUALITY_CHECKS_SUMMARY=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'report-summary.tsv')

rule FastQC_additional_checks:
    """
    Tests additional quality metrics based on the FastQC data file output.
    """
    input:
        TXT=TRIMMED_READS_QC_TXT,
        TXT_RAW=ORIG_READS_QC_TXT
    output:
        INFORMS=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'fastqc_checks', 'informs.io')
    params:
        running_dir=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'fastqc_checks')
    run:
        from camel.app.tools.fastqc.fastqcadditionalchecks import FastQCAdditionalChecks
        fastqc_checks = FastQCAdditionalChecks(camel)
        step = Step(rule, fastqc_checks, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(fastqc_checks, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastqc_checks, output)

rule Assembly_indexing:
    """
    Creates a Bowtie2 index for the assembly.
    """
    input:
        FASTA_REF=FASTA_ASSEMBLY
    output:
        INDEX_GENOME_PREFIX=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'assembly_indexing', 'genome_prefix.io')
    params:
        running_dir=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'assembly_indexing')
    run:
        from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index
        bowtie2_index = Bowtie2Index(camel)
        step = Step(rule, bowtie2_index, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(bowtie2_index, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bowtie2_index, output)

rule Read_mapping:
    """
    Maps the trimmed reads to the assembly.
    """
    input:
        FASTQ_PE=TRIMMED_READS_PE,
        FASTQ_SE_FORWARD=TRIMMED_READS_SE_FORWARD,
        FASTQ_SE_REVERSE=TRIMMED_READS_SE_REVERSE,
        INDEX_GENOME_PREFIX=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'assembly_indexing', 'genome_prefix.io')
    output:
        SAM=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'read_mapping', 'sam.io'),
        INFORMS=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'read_mapping', 'informs.io')
    params:
        running_dir=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'read_mapping')
    run:
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        bowtie2_map = Bowtie2Map(camel)
        step = Step(rule, bowtie2_map, camel, params.running_dir, config)
        # set proper maximum_fragment_length for 2x300 sequencing (MiSeq)
        bowtie2_map.update_parameters(maximum_fragment_length=1500)
        SnakemakeUtils.add_pickle_inputs(bowtie2_map, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bowtie2_map, output)

rule Sam_to_bam_conversion:
    """
    Converts the mapped reads SAM file to BAM format.
    """
    input:
        SAM=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'read_mapping', 'sam.io')
    output:
        BAM=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'read_mapping', 'bam.io')
    params:
        running_dir=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'read_mapping')
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        samtools_view = SamtoolsView(camel)
        step = Step(rule, samtools_view, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_view, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_view, output)

rule Alignment_sorting:
    """
    Sorts the alignment.
    """
    input:
        BAM=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'read_mapping', 'bam.io')
    output:
        BAM=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'alignment_sorting', 'bam-sorted.io')
    params:
        running_dir=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'alignment_sorting')
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort(camel)
        step = Step(rule, samtools_sort, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_sort, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_sort, output)

rule Coverage_calculation:
    """
    Determines the read depth at every position to estimate the coverage.
    """
    input:
        BAM=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'alignment_sorting', 'bam-sorted.io')
    output:
        INFORMS=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'coverage_calculation', 'informs.io')
    params:
        running_dir=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'coverage_calculation')
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth(camel)
        step = Step(rule, samtools_depth, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)

rule Get_cgMLST_stats:
    """
    Retrieves the number of cgMLST genes that were detected (only perfect hits are considered).
    """
    input:
        hits=os.path.join(TYPING_WORKING_DIR, 'cgMLST', 'hits-combined.io')
    output:
        INFORMS=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'informs-cgmlst.io')
    run:
        all_hits = SnakemakeUtils.load_object(input.hits)
        nb_perfect = len([v for v in all_hits if v.value.is_perfect_hit()])
        SnakemakeUtils.dump_object({'hits_found': nb_perfect, 'nb_of_loci': len(all_hits)}, output.INFORMS)

rule Get_cgMLST_dummy_stats:
    """
    Dummy rule to create dummy cgMLST stats when cgmlst is skipped
    """
    output:
        INFORMS=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'informs-cgmlst-dummy.io')
    run:
        SnakemakeUtils.dump_object({'hits_found': 1, 'nb_of_loci': 1}, output.INFORMS)
        
rule Check_quality_criteria:
    """
    Checks if the quality criteria were met.
    """
    input:
        INFORMS_cgmlst=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'informs-cgmlst.io') if 'cgMLST' in config['sequence_typing'] else os.path.join(QUALITY_CHECKS_WORKING_DIR, 'informs-cgmlst-dummy.io'),
        INFORMS_coverage=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'coverage_calculation', 'informs.io'),
        INFORMS_fastqc_checks=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'fastqc_checks', 'informs.io'),
        INFORMS_mapping=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'read_mapping', 'informs.io')
    output:
        INFORMS=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'quality_checks', 'informs.io')
    params:
        running_dir=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'quality_checks')
    run:
        from camel.app.tools.pipelines.quality_checks.qualitycriteriachecker import QualityCriteriaChecker
        quality_checker = QualityCriteriaChecker(camel)
        step = Step(rule, quality_checker, camel, params.running_dir, config)
        # cgMLST detection perc. cutoff: 90 minimal requirement (pandemic), 95 acceptable (isolates), 97 good
        quality_checker.update_parameters(minimal_cgmlst_percentage_genes_fail=90, minimal_cgmlst_percentage_genes_warn=95)
        SnakemakeUtils.add_pickle_inputs(quality_checker, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quality_checker, output)

rule Report_quality_checks:
    """
    Creates a report containing the quality checks information.
    """
    input:
        INFORMS_coverage=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'coverage_calculation', 'informs.io'),
        INFORMS_mapping=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'read_mapping', 'informs.io'),
        INFORMS_additional_checks=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'fastqc_checks', 'informs.io'),
        INFORMS_quality_criteria=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'quality_checks', 'informs.io'),
        INFORMS_cgmlst=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'informs-cgmlst.io') if 'cgMLST' in config['sequence_typing'] else os.path.join(QUALITY_CHECKS_WORKING_DIR, 'informs-cgmlst-dummy.io')
    output:
        VAL_HTML=QUALITY_CHECKS_REPORT
    params:
        running_dir=os.path.join(QUALITY_CHECKS_WORKING_DIR)
    run:
        from camel.app.tools.pipelines.quality_checks.htmlreporterqualitychecks import HtmlReporterQualityChecks
        reporter = HtmlReporterQualityChecks(camel)
        step = Step(rule, reporter, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule Summary_quality_checks:
    """
    Creates a tabular summary containing the quality checks information.
    """
    input:
        INFORMS_coverage=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'coverage_calculation', 'informs.io'),
        INFORMS_mapping=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'read_mapping', 'informs.io'),
        INFORMS_additional_checks=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'fastqc_checks', 'informs.io'),
        INFORMS_cgmlst=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'informs-cgmlst.io') if 'cgMLST' in config['sequence_typing'] else os.path.join(QUALITY_CHECKS_WORKING_DIR, 'informs-cgmlst-dummy.io')
    output:
        QUALITY_CHECKS_SUMMARY
    params:
        running_dir=os.path.join(QUALITY_CHECKS_WORKING_DIR)
    run:
        informs_cgmlst = SnakemakeUtils.load_object(input.INFORMS_cgmlst)
        informs_qc_checks = SnakemakeUtils.load_object(input.INFORMS_additional_checks)['tests']
        summary_data = [
            ('coverage', SnakemakeUtils.load_object(input.INFORMS_coverage)['median_depth']),
            ('assembly_mapping_rate', SnakemakeUtils.load_object(input.INFORMS_mapping)['stats_map_rate']),
            ('cgmlst_hits_found', informs_cgmlst['hits_found']),
            ('cgmlst_nb_loci', informs_cgmlst['nb_of_loci']),
            ('check_avg_qs_fwd', informs_qc_checks['Average quality score'][0]),
            ('check_avg_qs_rev', informs_qc_checks['Average quality score'][1]),
            ('check_gc_fwd', informs_qc_checks['GC content'][0]),
            ('check_gc_rev', informs_qc_checks['GC content'][1]),
            ('check_n_fraction_fwd', informs_qc_checks['Maximal N-fraction'][0]),
            ('check_n_fraction_rev', informs_qc_checks['Maximal N-fraction'][1]),
            ('check_qs_drop_fwd', informs_qc_checks['Mean Q-score drop'][0]),
            ('check_qs_drop_rev', informs_qc_checks['Mean Q-score drop'][1]),
            ('check_pb_seq_content_fwd', informs_qc_checks['Per base sequence content'][0]),
            ('check_pb_seq_content_rev', informs_qc_checks['Per base sequence content'][1]),
            ('check_slen_dist_fwd', informs_qc_checks['Sequence length distribution'][0]),
            ('check_slen_dist_rev', informs_qc_checks['Sequence length distribution'][1]),
        ]
        with open(output[0], 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')
