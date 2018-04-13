import os

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

camel = Camel()
working_dir = config['working_dir']

rule all:
    # This rule makes sure that all other rules are executed.
    input:
        os.path.join(working_dir, "report_read_trimming/html-report.io")

rule prepare_initial_input:
    # In this rule the Snakemake input files are converted to Camel IO pickles.
    # The FASTQ and FASTA files will be converted to ToolIOFile objects, the DIR_HTML to a ToolIODirectory. Values can
    # be converted by passing them as parameters (not shown here).
    # Doing all the conversions in the beginning makes the following rules easier, because you only have to work with
    # Camel IO pickles.
    # Make sure that the output keys match the input keys (otherwise the keys argument has to be set).
    input:
        FASTQ=config['fastq_pe'],
        DIR_HTML=[os.path.dirname(config['report'])]
    output:
        FASTQ=os.path.join(working_dir, "initial_input/fastq.io"),
        DIR_HTML=os.path.join(working_dir, "initial_input/dir_html.io")
    run:
        SnakemakeUtils.pickle_snake_input(input, output)

rule fastqc_pre_trimming:
    # This rule creates a HTML report to assess raw data quality before the read trimming step.
    # The add_pickle_input is used to convert and add the FASTQ input to the tool.
    # After the tool is executed, the 'TXT' and 'HTML' outputs are dumped into Camel IO pickles using the
    # dump_tool_output function.
    input:
        FASTQ=os.path.join(working_dir, "initial_input/fastq.io")
    output:
        TXT=os.path.join(working_dir, "fastqc_pre_trimming/txt.io"),
        HTML=os.path.join(working_dir, "fastqc_pre_trimming/html.io")
    params:
        working_dir=os.path.join(working_dir, "fastqc_pre_trimming")
    threads: 8
    run:
        from camel.app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC(camel)
        fastqc.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(fastqc, 'FASTQ', input.FASTQ)
        fastqc.run(params.working_dir)
        SnakemakeUtils.dump_tool_output(fastqc, 'TXT', output.TXT)
        SnakemakeUtils.dump_tool_output(fastqc, 'HTML', output.HTML)


rule read_trimming:
    # In this rule, the reads are trimmed using Trimmomatic.
    # The same techniques as in the last rule are used, but the add_pickle_inputs is used. This function has the
    # same purpose as the add_pickle_input, but it loops over all keys of the snake input.
    # The dump_tool_outputs does the same for the tool outputs. It also stores the informs in a Camel IO pickle.
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
        from camel.app.tools.trimmomatic.trimmomatic import Trimmomatic
        trimmomatic = Trimmomatic(camel)
        trimmomatic.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_inputs(trimmomatic, input)
        trimmomatic.run(params.working_dir)
        SnakemakeUtils.dump_tool_outputs(trimmomatic, output)

rule fastqc_post_trimming:
    # In this rule reports are generated for the trimmed reads.
    # The run_tool function is used, this function will first call the add_pickle_inputs before and the
    # dump_tool_outputs after tool execution. You can also specify the working directory of the tool. It is recommended
    # to use this function to run CAMEL tools.
    input:
        FASTQ=os.path.join(working_dir, "read_trimming/fastq_pe.io")
    output:
        TXT=os.path.join(working_dir, "fastqc_post_trimming/txt.io"),
        HTML=os.path.join(working_dir, "fastqc_post_trimming/html.io")
    params:
        working_dir=os.path.join(working_dir, "fastqc_post_trimming")
    threads: 8
    run:
        from camel.app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC(camel)
        fastqc.update_parameters(threads=threads)
        SnakemakeUtils.run_tool(fastqc, input, output, params.working_dir)

rule report_read_trimming:
    # This rule creates the HTML report for the trimming pipeline.
    # The HTML input cannot be added as a pickle because it is a file that has to be created on the fly. The regular
    # Camel way of adding tool inputs is used instead. The other files are added the same way as before.
    # We also add the informs generated in the read trimming step. Keys that start with 'INFORMS' will be added as input
    # informs. In this case the informs are added with key 'trimming'.
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
        from camel.app.tools.pipelines.read_trimming.htmlreporterreadtrimming import HtmlReporterReadTrimming
        reporter = HtmlReporterReadTrimming(camel)
        report_path = os.path.join(params.working_dir, 'report_trimming.html')
        open(report_path, 'w').close()
        reporter.add_input_files({'HTML': [ToolIOFile(report_path)]})
        SnakemakeUtils.run_tool(reporter, input, output, params.working_dir)
