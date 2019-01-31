import os

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.read_trimming_iontorrent import OUTPUT_TRIMMING_IT_READS, OUTPUT_TRIMMING_IT_REPORT, \
    OUTPUT_TRIMMING_IT_SUMMARY

camel = Camel.get_instance()

rule Timming_iontorrent_pickle_fastq_input:
    """
    Creates a pickle for the fastq input files.
    """
    input:
        FASTQ=config['fastq_se'][0]['path'] if config.get('read_type', 'illumina') == 'iontorrent' else []
    output:
        FASTQ_SE=os.path.join(config['working_dir'], 'read_trimming_it', 'input', 'fastq-se.io')
    run:
        from camel.app.io.tooliofile import ToolIOFile
        SnakemakeUtils.dump_object([ToolIOFile(input.FASTQ)], output.FASTQ_SE)

rule Trimming_iontorrent_fastqc_pre:
    """
    Creates FastQC reports for the raw reads. 
    """
    input:
        FASTQ=os.path.join(config['working_dir'], 'read_trimming_it', 'input', 'fastq-se.io')
    output:
        HTML=os.path.join(config['working_dir'], 'read_trimming_it', 'fastqc-pre', 'html.io'),
        TXT=os.path.join(config['working_dir'], 'read_trimming_it', 'fastqc-pre', 'txt.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'read_trimming_it', 'fastqc-pre')
    threads: 2
    run:
        from camel.app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC(camel)
        SnakemakeUtils.add_pickle_inputs(fastqc, input)
        step = Step(rule, fastqc, camel, params.running_dir, config)
        fastqc.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastqc, output)

rule Trimming_iontorrent_filter_length:
    """
    Filters input reads based on read length.
    """
    input:
        FASTQ=os.path.join(config['working_dir'], 'read_trimming_it', 'input', 'fastq-se.io')
    output:
        FASTQ=os.path.join(config['working_dir'], 'read_trimming_it', 'trim_length', 'fastq.io'),
        INFORMS=os.path.join(config['working_dir'], 'read_trimming_it', 'trim_length', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'read_trimming_it', 'trim_length')
    run:
        from camel.app.tools.fastx.fastqqualitytrimmer import FastqQualityTrimmer
        trimmer = FastqQualityTrimmer(camel)
        SnakemakeUtils.add_pickle_inputs(trimmer, input)
        step = Step(rule, trimmer, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(trimmer, output)

rule Trimming_iontorrent_filter_quality:
    """
    Filters input reads based on quality score.
    """
    input:
        FASTQ=os.path.join(config['working_dir'], 'read_trimming_it', 'trim_length', 'fastq.io')
    output:
        FASTQ=os.path.join(config['working_dir'], OUTPUT_TRIMMING_IT_READS),
        INFORMS=os.path.join(config['working_dir'], 'read_trimming_it', 'trim_qual', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'read_trimming_it', 'trim_qual')
    run:
        from camel.app.tools.fastx.fastqqualityfilter import FastqQualityFilter
        q_filter = FastqQualityFilter(camel)
        SnakemakeUtils.add_pickle_inputs(q_filter , input)
        step = Step(rule, q_filter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(q_filter, output)

rule Trimming_iontorrent_fastqc_post:
    """
    Creates FastQC reports of the filtered reads.
    """
    input:
        FASTQ=os.path.join(config['working_dir'], 'read_trimming_it', 'trim_qual', 'fastq.io')
    output:
        HTML=os.path.join(config['working_dir'], 'read_trimming_it', 'fastqc-post', 'html.io'),
        TXT=os.path.join(config['working_dir'], 'read_trimming_it', 'fastqc-post', 'txt.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'read_trimming_it', 'fastqc-post')
    threads: 4
    run:
        from camel.app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC(camel)
        SnakemakeUtils.add_pickle_inputs(fastqc, input)
        step = Step(rule, fastqc, camel, params.running_dir, config)
        fastqc.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastqc, output)

rule Trimming_iontorrent_report:
    """
    Creates the HTML report with the trimming output files and statistics.
    """
    input:
        HTML_Pre=os.path.join(config['working_dir'], 'read_trimming_it', 'fastqc-pre', 'html.io'),
        HTML_Post=os.path.join(config['working_dir'], 'read_trimming_it', 'fastqc-post', 'html.io'),
        FASTQ=os.path.join(config['working_dir'], OUTPUT_TRIMMING_IT_READS),
        INFORMS_filt_len=os.path.join(config['working_dir'], 'read_trimming_it', 'trim_length', 'informs.io'),
        INFORMS_filt_qual=os.path.join(config['working_dir'], 'read_trimming_it', 'trim_qual', 'informs.io')
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_TRIMMING_IT_REPORT)
    params:
        running_dir=os.path.join(config['working_dir'], 'read_trimming_it', 'report')
    run:
        from camel.app.tools.pipelines.read_trimming.reportertrimmingiontorrent import ReporterTrimmingIonTorrent
        reporter = ReporterTrimmingIonTorrent(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule Trimming_iontorrent_collect_summary_info:
    """
    Collects the summary information for the IonTorrent read trimming.
    """
    input:
        INFORMS_filt_len=os.path.join(config['working_dir'], 'read_trimming_it', 'trim_length', 'informs.io'),
        INFORMS_filt_qual=os.path.join(config['working_dir'], 'read_trimming_it', 'trim_qual', 'informs.io')
    output:
        os.path.join(config['working_dir'], OUTPUT_TRIMMING_IT_SUMMARY)
    run:
        informs_len = SnakemakeUtils.load_object(input.INFORMS_filt_len)
        informs_qual = SnakemakeUtils.load_object(input.INFORMS_filt_qual)
        with open(output[0], 'w') as handle:
            for k, v in [
                ('filt_len_in', informs_len['input_reads']),
                ('filt_len_out', informs_len['output_reads']),
                ('filt_qual_in', informs_qual['input_reads']),
                ('filt_qual_out', informs_qual['output_reads'])]:
                handle.write('\t'.join([k, str(v)]))
                handle.write('\n')
