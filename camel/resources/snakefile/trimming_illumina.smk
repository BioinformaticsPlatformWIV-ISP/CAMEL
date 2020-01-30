from pathlib import Path

from camel.app.camel import Camel
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import trimming_illumina
from camel.app.pipeline.step import Step

camel = Camel.get_instance()


rule trimming_illumina_pickle_input:
    """
    Creates pickled FASTQ PE input.
    """
    input:
        FASTQ = [x['path'] for x in config.get('fastq_pe', [])]
    output:
        FASTQ_PE = Path(config['working_dir']) / 'trimming_illumina' / 'input'/ 'fastq-pe.io'
    run:
        from camel.app.io.tooliofile import ToolIOFile
        SnakemakeUtils.dump_object([ToolIOFile(x) for x in input.FASTQ], output.FASTQ_PE)


rule trimming_illumina_fastqc_pre:
    """
    Creates FastQC reports for the raw reads. 
    """
    input:
        FASTQ = rules.trimming_illumina_pickle_input.output.FASTQ_PE
    output:
        HTML = Path(config['working_dir']) / 'trimming_illumina' / 'fastqc-pre' / 'html.io',
        TXT = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_FASTQC_PRE
    params:
        running_dir = Path(config['working_dir']) / 'trimming_illumina' / 'fastqc-pre'
    threads: 4
    run:
        from camel.app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC(camel)
        SnakemakeUtils.add_pickle_inputs(fastqc, input)
        step = Step(rule, fastqc, camel, params.running_dir, config)
        fastqc.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastqc, output)


rule trimming_illumina_trimmomatic:
    """
    Read trimming using trimmomatic.
    """
    input:
        FASTQ_PE = rules.trimming_illumina_pickle_input.output.FASTQ_PE
    output:
        FASTQ_PE = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_READS_PE,
        FASTQ_SE_FORWARD = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_READS_SE_FWD,
        FASTQ_SE_REVERSE = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_READS_SE_REV,
        INFORMS = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_INFORMS
    threads: 4
    priority: 1
    params:
        running_dir = Path(config['working_dir']) / 'trimming_illumina' / 'trimmomatic'
    run:
        from camel.app.tools.trimmomatic.trimmomatic import Trimmomatic
        trimmomatic = Trimmomatic(camel)
        SnakemakeUtils.add_pickle_inputs(trimmomatic, input)
        step = Step(rule, trimmomatic, camel, params.running_dir, config)
        trimmomatic.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(trimmomatic, output)


rule trimming_illumina_fastqc_post:
    """
    Creates FastQC reports of the trimmed reads.
    """
    input:
        FASTQ = rules.trimming_illumina_trimmomatic.output.FASTQ_PE
    output:
        HTML = Path(config['working_dir']) /  'trimming_illumina' / 'fastqc-post' / 'html.io',
        TXT = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_FASTQC_POST
    params:
        running_dir = Path(config['working_dir']) / 'trimming_illumina' / 'fastqc-post'
    threads: 4
    run:
        from camel.app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC(camel)
        SnakemakeUtils.add_pickle_inputs(fastqc, input)
        step = Step(rule, fastqc, camel, params.running_dir, config)
        fastqc.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastqc, output)


rule trimming_illumina_report:
    """
    Creates the HTML report with the trimming output files and statistics.
    """
    input:
        HTML_PRE = rules.trimming_illumina_fastqc_pre.output.HTML,
        HTML_POST = rules.trimming_illumina_fastqc_post.output.HTML,
        FASTQ_PE = rules.trimming_illumina_trimmomatic.output.FASTQ_PE,
        FASTQ_SE_FORWARD = rules.trimming_illumina_trimmomatic.output.FASTQ_SE_FORWARD,
        FASTQ_SE_REVERSE = rules.trimming_illumina_trimmomatic.output.FASTQ_SE_REVERSE,
        INFORMS_trimming = rules.trimming_illumina_trimmomatic.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'trimming_illumina' / 'report',
        export_fastq = config['read_trimming'].get('export_fastq', 'false')
    run:
        from camel.app.tools.pipelines.read_trimming.reportertrimming import ReporterTrimming
        reporter = ReporterTrimming(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.running_dir, config)
        reporter.update_parameters(export_fastq=params.export_fastq)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)


rule trimming_illumina_dump_summary_info:
    """
    Dumps the summary information from the read trimming pipeline.
    """
    input:
        INFORMS_trimming = rules.trimming_illumina_trimmomatic.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_SUMMARY
    params:
        running_dir = Path(config['working_dir']) / 'trimming_illumina' / 'summary'
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


rule trimming_illumina_to_dict:
    """
    Combines the trimmed reads into a dictionary.
    """
    input:
        FASTQ_PE = rules.trimming_illumina_trimmomatic.output.FASTQ_PE,
        FASTQ_SE_FWD = rules.trimming_illumina_trimmomatic.output.FASTQ_SE_FORWARD,
        FASTQ_SE_REV = rules.trimming_illumina_trimmomatic.output.FASTQ_SE_REVERSE
    output:
        IO = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_DICT
    run:
        output_dict = {
            'PE': SnakemakeUtils.load_object(input.FASTQ_PE)
        }
        se_fwd = SnakemakeUtils.load_object(input.FASTQ_SE_FWD)
        if len(se_fwd) > 0:
            output_dict['SE_FWD'] = se_fwd
        se_rev = SnakemakeUtils.load_object(input.FASTQ_SE_REV)
        if len(se_rev) > 0:
            output_dict['SE_REV'] = se_rev
        SnakemakeUtils.dump_object(output_dict, output.IO)
