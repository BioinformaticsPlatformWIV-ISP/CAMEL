rule fastqc_pre_trimming:
    """
    Creates FastQC reports for the raw reads.
    """
    input:
        FASTQ = config['fastq_pe']
    output:
        HTML = os.path.join(__WORKING_DIR, 'pre_trimming', 'html.io'),
        TXT = os.path.join(__WORKING_DIR, 'pre_trimming', 'txt.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'pre_trimming')
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
        FASTQ_PE = os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-pe.io'),
        FASTQ_SE_FORWARD = os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-se-forward.io'),
        FASTQ_SE_REVERSE = os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-se-reverse.io'),
        INFORMS = os.path.join(__WORKING_DIR, 'read_trimming', 'informs.io')
    threads:
        8
    params:
        running_dir = os.path.join(__WORKING_DIR, 'read_trimming')
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
        FASTQ = os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-pe.io'),
    output:
        HTML = os.path.join(__WORKING_DIR, 'post_trimming', 'html.io'),
        TXT = os.path.join(__WORKING_DIR, 'post_trimming', 'txt.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'post_trimming')
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
        HTML_Pre = os.path.join(__WORKING_DIR, 'pre_trimming', 'html.io'),
        HTML_Post = os.path.join(__WORKING_DIR, 'post_trimming', 'html.io'),
        FASTQ_PE = os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-pe.io'),
        FASTQ_SE_FORWARD = os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-se-forward.io'),
        FASTQ_SE_REVERSE = os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-se-reverse.io'),
        INFORMS_trimming = os.path.join(__WORKING_DIR, 'read_trimming', 'informs.io')
    output:
        VAL_HTML = os.path.join(__WORKING_DIR, 'report_read_trimming', 'html.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'report_read_trimming'),
        output_dir = config['output_dir']
    run:
        from app.tools.pipelines.read_trimming.htmlreporterreadtrimming import HtmlReporterReadTrimming
        reporter = HtmlReporterReadTrimming(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input, optionals=['FASTQ_SE_FORWARD', 'FASTQ_SE_REVERSE'])
        step = SnakeStep(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        reporter.tool_outputs['VAL_HTML'][0].value.copy_files(params.output_dir)
        SnakemakeUtils.dump_tool_outputs(reporter, output)
