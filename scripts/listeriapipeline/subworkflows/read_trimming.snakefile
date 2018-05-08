READ_TRIMMING_WORKING_DIR = os.path.join(__WORKING_DIR, 'read_trimming')
READ_TRIMMING_REPORT = os.path.join(READ_TRIMMING_WORKING_DIR, 'report-html.io')
READ_TRIMMING_SUMMARY = os.path.join(READ_TRIMMING_WORKING_DIR, 'report-summary.tsv')
TRIMMED_READS_PE = os.path.join(READ_TRIMMING_WORKING_DIR, 'trimming', 'fastq-pe.io')
TRIMMED_READS_SE_FORWARD = os.path.join(READ_TRIMMING_WORKING_DIR, 'trimming', 'fastq-se-forward.io')
TRIMMED_READS_SE_REVERSE = os.path.join(READ_TRIMMING_WORKING_DIR, 'trimming', 'fastq-se-reverse.io')
TRIMMING_INFORM = os.path.join(READ_TRIMMING_WORKING_DIR, 'trimming', 'informs.io')
ORIG_READS_QC_TXT = os.path.join(READ_TRIMMING_WORKING_DIR, 'pre_trimming', 'txt.io')
TRIMMED_READS_QC_TXT = os.path.join(READ_TRIMMING_WORKING_DIR, 'post_trimming', 'txt.io')

rule fastqc_pre_trimming:
    """
    Creates FastQC reports for the raw reads.
    """
    input:
        FASTQ = config['fastq_pe']
    output:
        HTML = os.path.join(READ_TRIMMING_WORKING_DIR, 'pre_trimming', 'html.io'),
        TXT = ORIG_READS_QC_TXT
    params:
        running_dir = os.path.join(READ_TRIMMING_WORKING_DIR, 'pre_trimming')
    threads:
        8
    run:
        from app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC(camel)
        fastqc.add_input_files({'FASTQ': [ToolIOFile(x) for x in input.FASTQ]})
        step = SnakeStep(rule, fastqc, camel, params.running_dir, config)
        fastqc.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastqc, output)

rule read_trimming:
    """
    Read trimming using trimmomatic.
    """
    input:
        FASTQ = config['fastq_pe']
    output:
        FASTQ_PE = TRIMMED_READS_PE,
        FASTQ_SE_FORWARD = TRIMMED_READS_SE_FORWARD,
        FASTQ_SE_REVERSE = TRIMMED_READS_SE_REVERSE,
        INFORMS = TRIMMING_INFORM
    threads:
        8
    params:
        running_dir = os.path.join(READ_TRIMMING_WORKING_DIR, 'trimming')
    run:
        from app.tools.trimmomatic.trimmomatic import Trimmomatic
        trimmomatic = Trimmomatic(camel)
        trimmomatic.add_input_files({'FASTQ_PE': [ToolIOFile(x) for x in input.FASTQ]})
        step = SnakeStep(rule, trimmomatic, camel, params.running_dir, config)
        trimmomatic.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(trimmomatic, output)

rule fastqc_post_trimming:
    """
    Creates FastQC reports of the trimmed reads.
    """
    input:
        FASTQ = TRIMMED_READS_PE
    output:
        HTML = os.path.join(READ_TRIMMING_WORKING_DIR, 'post_trimming', 'html.io'),
        TXT = TRIMMED_READS_QC_TXT
    params:
        running_dir = os.path.join(READ_TRIMMING_WORKING_DIR, 'post_trimming')
    threads:
        8
    run:
        from app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC(camel)
        SnakemakeUtils.add_pickle_inputs(fastqc, input)
        step = SnakeStep(rule, fastqc, camel, params.running_dir, config)
        fastqc.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastqc, output)

rule report_read_trimming:
    """
    Creates the HTML report with the trimming output files and statistics.
    """
    input:
        HTML_Pre = os.path.join(READ_TRIMMING_WORKING_DIR, 'pre_trimming', 'html.io'),
        HTML_Post = os.path.join(READ_TRIMMING_WORKING_DIR, 'post_trimming', 'html.io'),
        FASTQ_PE = TRIMMED_READS_PE,
        FASTQ_SE_FORWARD = TRIMMED_READS_SE_FORWARD,
        FASTQ_SE_REVERSE = TRIMMED_READS_SE_REVERSE,
        INFORMS_trimming = TRIMMING_INFORM
    output:
        VAL_HTML = READ_TRIMMING_REPORT
    params:
        running_dir = os.path.join(READ_TRIMMING_WORKING_DIR),
        output_dir = config['output_dir']
    run:
        from app.tools.pipelines.read_trimming.htmlreporterreadtrimming import HtmlReporterReadTrimming
        reporter = HtmlReporterReadTrimming(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input, optionals=['FASTQ_SE_FORWARD', 'FASTQ_SE_REVERSE'])
        step = SnakeStep(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        reporter.tool_outputs['VAL_HTML'][0].value.copy_files(params.output_dir)
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule summary_read_trimming:
    """
    Creates a tabular summary with the trimming output files and statistics.
    """
    input:
        INFORMS_trimming=os.path.join(READ_TRIMMING_WORKING_DIR, 'trimming', 'informs.io')
    output:
        READ_TRIMMING_SUMMARY
    params:
        running_dir=os.path.join(config['working_dir'], 'read_trimming', 'summary')
    run:
        trimmomatic_informs = SnakemakeUtils.load_object(input.INFORMS_trimming)
        summary_data = [
            ('total_reads_pairs', trimmomatic_informs['paired_reads_in']),
            ('total_reads_pairs_trimmed', trimmomatic_informs['paired_reads_out'].split(' ')[0]),
            ('forward_only_surviving', trimmomatic_informs['forward_only_reads'].split(' ')[0]),
            ('reverse_only_surviving', trimmomatic_informs['reverse_only_reads'].split(' ')[0]),
            ('dropped', trimmomatic_informs['reads_drop'].split(' ')[0])
        ]
        with open(output[0], 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')
