import os

from app.camel import Camel
from app.io.tooliofile import ToolIOFile
from app.pipeline.snakestep import SnakeStep
from app.snakemake.snakemakeutils import SnakemakeUtils

camel = Camel()
working_dir = config['working_dir']

rule all:
    # This rule makes sure that all other rules are executed.
    input:
        os.path.join(working_dir, "report_read_trimming/html-report.io")

rule Prepare_initial_input:
    # This rule is the same as the other example, as this is not a part of the pipeline.
    input:
        FASTQ=config['fastq_pe'],
        DIR_HTML=[os.path.dirname(config['report'])]
    output:
        FASTQ=os.path.join(working_dir, "initial_input/fastq.io"),
        DIR_HTML=os.path.join(working_dir, "initial_input/dir_html.io")
    run:
        SnakemakeUtils.pickle_snake_input(input, output)

rule FastQC_pre_trimming:
    # This rule creates a HTML report to assess raw data quality before the read trimming step.
    # The same methods to link input / output are used but instead of running the tool directly, it is ran trough a
    # step object. This step object ensures that:
    # - Parameters are retrieved from the database
    # - Input / outputs are logged (IF it is enabled in the config)
    #
    # It is important that the rule names correspond to the step names in the database / YAML file.
    input:
        FASTQ=os.path.join(working_dir, "initial_input/fastq.io")
    output:
        TXT=os.path.join(working_dir, "fastqc_pre_trimming/txt.io"),
        HTML=os.path.join(working_dir, "fastqc_pre_trimming/html.io")
    params:
        working_dir=os.path.join(working_dir, "fastqc_pre_trimming")
    threads: 8
    run:
        from app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC(camel)
        SnakemakeUtils.add_pickle_inputs(fastqc, input)
        step = SnakeStep(rule, fastqc, camel, params.working_dir, config)
        #
        # update_parameters should happen after creating SnakeStep which retrieves the default parameters from DB
        fastqc.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastqc, output)

rule Read_trimming:
    # In this rule, the reads are trimmed using Trimmomatic.
    input:
        FASTQ_PE=os.path.join(working_dir, "initial_input/fastq.io")
    output:
        FASTQ_PE=os.path.join(working_dir, "read_trimming/fastq_pe.io"),
        FASTQ_SE_FORWARD=os.path.join(working_dir, "read_trimming/fastq_se_forward.io"),
        FASTQ_SE_REVERSE=os.path.join(working_dir, "read_trimming/fastq_se_reverse.io"),
        INFORMS=os.path.join(working_dir, "read_trimming/informs.io")
    params:
        working_dir=os.path.join(working_dir, "read_trimming")
    threads: 8
    run:
        from app.tools.trimmomatic.trimmomatic import Trimmomatic
        trimmomatic = Trimmomatic(camel)
        SnakemakeUtils.add_pickle_inputs(trimmomatic, input)
        step = SnakeStep(rule, trimmomatic, camel, params.working_dir, config)
        trimmomatic.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(trimmomatic, output)

rule FastQC_post_trimming:
    # In this rule quality reports are generated for the trimmed reads.
    input:
        FASTQ=os.path.join(working_dir, "read_trimming/fastq_pe.io")
    output:
        TXT=os.path.join(working_dir, "fastqc_post_trimming/txt.io"),
        HTML=os.path.join(working_dir, "fastqc_post_trimming/html.io")
    params:
        working_dir=os.path.join(working_dir, "fastqc_post_trimming")
    threads: 8
    run:
        from app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC(camel)
        SnakemakeUtils.add_pickle_inputs(fastqc, input)
        step = SnakeStep(rule, fastqc, camel, params.working_dir, config)
        fastqc.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastqc, output)

rule Report_generation:
    # This rule creates the HTML report for the trimming pipeline.
    input:
        HTML_Pre=os.path.join(working_dir, "fastqc_pre_trimming/html.io"),
        HTML_Post=os.path.join(working_dir, "fastqc_post_trimming/html.io"),
        FASTQ_PE=os.path.join(working_dir, "read_trimming/fastq_pe.io"),
        FASTQ_SE_FORWARD=os.path.join(working_dir, "read_trimming/fastq_se_forward.io"),
        FASTQ_SE_REVERSE=os.path.join(working_dir, "read_trimming/fastq_se_reverse.io"),
        DIR_HTML=os.path.join(working_dir, "initial_input/dir_html.io"),
        INFORMS_trimming=os.path.join(working_dir, "read_trimming/informs.io")
    output:
        HTML=os.path.join(working_dir, "report_read_trimming/html-report.io")
    params:
        working_dir = os.path.join(working_dir, "report_read_trimming")
    run:
        from app.tools.pipelines.read_trimming.htmlreporterreadtrimming import HtmlReporterReadTrimming
        reporter = HtmlReporterReadTrimming(camel)
        report_path = os.path.join(params.working_dir, 'report_trimming.html')
        open(report_path, 'w').close()
        reporter.add_input_files({'HTML': [ToolIOFile(report_path)]})
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = SnakeStep(rule, reporter, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.run_tool(reporter, input, output, params.working_dir)
        SnakemakeUtils.dump_tool_outputs(reporter, output)
