"""
This Snakefile is used for read trimming. It also generates a report with some trimming statistics and FASTQ reports
before and after trimming.
"""
import os

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_READS_PE, OUTPUT_READ_TRIMMING_READS_SE_FWD, \
    OUTPUT_READ_TRIMMING_READS_SE_REV, OUTPUT_READ_TRIMMING_REPORT, OUTPUT_READ_TRIMMING_SUMMARY, \
    OUTPUT_READ_TRIMMING_INFORMS, OUTPUT_READ_TRIMMING_FASTQC_POST, OUTPUT_READ_TRIMMING_FASTQC_PRE

camel = Camel.get_instance()


rule Pickle_fastq_input:
    """
    Creates pickled FASTQ PE input.
    """
    input:
        FASTQ=[x['path'] for x in config.get('fastq_pe', [])]
    output:
        FASTQ_PE=os.path.join(config['working_dir'], 'read_trimming', 'input', 'fastq-pe.io')
    run:
        SnakemakeUtils.dump_object([ToolIOFile(x) for x in input.FASTQ], output.FASTQ_PE)

rule Trimming_fastqc_pre:
    """
    Creates FastQC reports for the raw reads. 
    """
    input:
        FASTQ=os.path.join(config['working_dir'], 'read_trimming', 'input', 'fastq-pe.io')
    output:
        HTML=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_FASTQC_PRE),
        TXT=os.path.join(config['working_dir'], 'read_trimming', 'fastqc-pre', 'txt.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'read_trimming', 'fastqc-pre')
    threads: 4
    run:
        from camel.app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC(camel)
        SnakemakeUtils.add_pickle_inputs(fastqc, input)
        step = Step(rule, fastqc, camel, params.running_dir, config)
        fastqc.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastqc, output)

rule Trimming_trimmomatic:
    """
    Read trimming using trimmomatic.
    """
    input:
        FASTQ_PE=os.path.join(config['working_dir'], 'read_trimming', 'input', 'fastq-pe.io')
    output:
        FASTQ_PE=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_PE),
        FASTQ_SE_FORWARD=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_SE_FWD),
        FASTQ_SE_REVERSE=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_SE_REV),
        INFORMS=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_INFORMS)
    threads: 4
    priority: 1
    params:
        running_dir=os.path.join(config['working_dir'], 'read_trimming', 'trimmomatic')
    run:
        from camel.app.tools.trimmomatic.trimmomatic import Trimmomatic
        trimmomatic = Trimmomatic(camel)
        SnakemakeUtils.add_pickle_inputs(trimmomatic, input)
        step = Step(rule, trimmomatic, camel, params.running_dir, config)
        trimmomatic.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(trimmomatic, output)

rule Trimming_fastqc_post:
    """
    Creates FastQC reports of the trimmed reads.
    """
    input:
        FASTQ=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_PE)
    output:
        HTML=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_FASTQC_POST),
        TXT=os.path.join(config['working_dir'], 'read_trimming', 'fastqc-post', 'txt.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'read_trimming', 'fastqc-post')
    threads: 4
    run:
        from camel.app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC(camel)
        SnakemakeUtils.add_pickle_inputs(fastqc, input)
        step = Step(rule, fastqc, camel, params.running_dir, config)
        fastqc.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastqc, output)

rule Trimming_report:
    """
    Creates the HTML report with the trimming output files and statistics.
    """
    input:
        HTML_Pre=os.path.join(config['working_dir'], 'read_trimming', 'fastqc-pre', 'html.io'),
        HTML_Post=os.path.join(config['working_dir'], 'read_trimming', 'fastqc-post', 'html.io'),
        FASTQ_PE=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_PE),
        FASTQ_SE_FORWARD=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_SE_FWD),
        FASTQ_SE_REVERSE=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_SE_REV),
        INFORMS_trimming=os.path.join(config['working_dir'], 'read_trimming', 'trimmomatic', 'informs.io')
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_REPORT)
    params:
        running_dir=os.path.join(config['working_dir'], 'read_trimming', 'report'),
        export_fastq=config['read_trimming'].get('export_fastq', 'false')
    run:
        from camel.app.tools.pipelines.read_trimming.reportertrimming import ReporterTrimming
        reporter = ReporterTrimming(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.running_dir, config)
        reporter.update_parameters(export_fastq=params.export_fastq)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule Trimming_dump_summary_info:
    """
    Dumps the summary information from the read trimming pipeline.
    """
    input:
        INFORMS_trimming=os.path.join(config['working_dir'], 'read_trimming', 'trimmomatic', 'informs.io')
    output:
        os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_SUMMARY)
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
