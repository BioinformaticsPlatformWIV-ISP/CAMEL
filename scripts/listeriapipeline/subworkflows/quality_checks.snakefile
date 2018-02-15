QUALITY_CHECKS_WORKING_DIR=os.path.join(__WORKING_DIR, 'quality_checks')
QUALITY_CHECKS_REPORT=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'report-html.io')


rule FastQC_additional_checks:
    """
    Tests additional quality metrics based on the FastQC data file output.
    """
    input:
        TXT=TRIMMED_READS_QC_TXT
    output:
        INFORMS=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'fastqc_checks', 'informs.io')
    params:
        running_dir=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'fastqc_checks')
    run:
        from app.tools.fastqc.fastqcadditionalchecks import FastQCAdditionalChecks
        fastqc_checks = FastQCAdditionalChecks(camel)
        step = SnakeStep(rule, fastqc_checks, camel, params.running_dir, config)
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
        from app.tools.bowtie2.bowtie2index import Bowtie2Index
        bowtie2_index = Bowtie2Index(camel)
        step = SnakeStep(rule, bowtie2_index, camel, params.running_dir, config)
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
        from app.tools.bowtie2.bowtie2map import Bowtie2Map
        bowtie2_map = Bowtie2Map(camel)
        step = SnakeStep(rule, bowtie2_map, camel, params.running_dir, config)
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
        from app.tools.samtools.samtoolsview import SamtoolsView
        samtools_view = SamtoolsView(camel)
        step = SnakeStep(rule, samtools_view, camel, params.running_dir, config)
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
        from app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort(camel)
        step = SnakeStep(rule, samtools_sort, camel, params.running_dir, config)
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
        from app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth(camel)
        step = SnakeStep(rule, samtools_depth, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)

rule Check_quality_criteria:
    """
    Checks if the quality criteria were met.
    """
    input:
        INFORMS_coverage=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'coverage_calculation', 'informs.io'),
        INFORMS_fastqc_checks=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'fastqc_checks', 'informs.io'),
        INFORMS_mapping=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'read_mapping', 'informs.io')
    output:
        INFORMS=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'quality_checks', 'informs.io')
    params:
        running_dir=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'quality_checks')
    run:
        from app.tools.pipelines.quality_checks.qualitycriteriachecker import QualityCriteriaChecker
        quality_checker = QualityCriteriaChecker(camel)
        step = SnakeStep(rule, quality_checker, camel, params.running_dir, config)
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
        INFORMS_quality_criteria=os.path.join(QUALITY_CHECKS_WORKING_DIR, 'quality_checks', 'informs.io')
    output:
        VAL_HTML=QUALITY_CHECKS_REPORT
    params:
        running_dir=os.path.join(QUALITY_CHECKS_WORKING_DIR)
    run:
        from app.tools.pipelines.quality_checks.htmlreporterqualitychecks import HtmlReporterQualityChecks
        reporter = HtmlReporterQualityChecks(camel)
        step = SnakeStep(rule, reporter, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)
