from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import trimming_ont

camel = Camel.get_instance()


rule trimming_ont_nanoplot_pre:
    """
    Creates NanoPlot reports for the raw reads. 
    """
    input:
        FASTQ = Path(config['working_dir']) / trimming_ont.INPUT_ONT_FASTQ
    output:
        HTML = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_HTML_PRE,
        TSV = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_TXT_PRE,
        INFORMS = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_INFORMS_PRE
    params:
        running_dir = Path(config['working_dir']) / 'trimming_ont' / 'nanoplot-pre'
    threads: 4
    run:
        from camel.app.tools.nanoplot.nanoplot import NanoPlot
        nanoplot = NanoPlot(camel)
        SnakemakeUtils.add_pickle_inputs(nanoplot, input)
        step = Step(str(rule), nanoplot, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(nanoplot, output)

rule trimming_ont_seqkit:
    """
    Read trimming using seqkit.
    """
    input:
        FASTQ = Path(config['working_dir']) / trimming_ont.INPUT_ONT_FASTQ
    output:
        FASTQ = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_READS,
        INFORMS = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_INFORMS
    threads: 4
    priority: 1
    params:
        running_dir = Path(config['working_dir']) / 'trimming_ont' / 'seqkit',
        sample_name = config.get('sample_name', 'reads'),
        min_length = config.get('read_trimming', {}).get('ont', {}).get('min_length', 500),
        min_qual = config.get('read_trimming', {}).get('ont', {}).get('min_qual', 7),
    run:
        from camel.app.tools.seqkit.seqkitseq import SeqkitSeq
        seqkit = SeqkitSeq(Camel.get_instance())
        seqkit.update_parameters(
            min_length=params.min_length,
            min_qual=params.min_qual,
            output_filename=f'{params.sample_name}-filtered.fastq',
            threads=threads
        )
        SnakemakeUtils.add_pickle_inputs(seqkit, input)
        step = Step(str(rule), seqkit, Camel.get_instance(), params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqkit, output)

rule trimming_ont_nanoplot_post:
    """
    Creates NanoPlot reports of the trimmed reads.
    """
    input:
        FASTQ = rules.trimming_ont_seqkit.output.FASTQ
    output:
        HTML = Path(config['working_dir']) /  trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_HTML_POST,
        TSV = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_TXT_POST,
        INFORMS = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_INFORMS_POST
    params:
        running_dir = Path(config['working_dir']) / 'trimming_ont' / 'nanoplot-post'
    threads: 4
    run:
        from camel.app.tools.nanoplot.nanoplot import NanoPlot
        nanoplot = NanoPlot(camel)
        SnakemakeUtils.add_pickle_inputs(nanoplot, input)
        step = Step(str(rule), nanoplot, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(nanoplot, output)

rule trimming_ont_report:
    """
    Creates the HTML report with the trimming output files and statistics.
    """
    input:
        HTML_PRE = rules.trimming_ont_nanoplot_pre.output.HTML,
        HTML_POST = rules.trimming_ont_nanoplot_post.output.HTML,
        FASTQ_PE = rules.trimming_ont_seqkit.output.FASTQ,
        INFORMS_trimming = rules.trimming_ont_seqkit.output.INFORMS,
        INFORMS_nanoplot_pre = rules.trimming_ont_nanoplot_pre.output.INFORMS,
        INFORMS_nanoplot_post= rules.trimming_ont_nanoplot_post.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'trimming_ont' / 'report',
        export_fastq = config['read_trimming'].get('export_fastq', 'false')
    run:
        from camel.app.tools.pipelines.read_trimming.reportertrimmingont import ReporterTrimmingONT
        reporter = ReporterTrimmingONT(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, camel, params.running_dir)
        reporter.update_parameters(export_fastq=str(params.export_fastq))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule trimming_ont_dump_summary_info:
    """
    Dumps the summary information from the read trimming pipeline.
    """
    input:
        INFORMS_trimming = rules.trimming_ont_seqkit.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_SUMMARY
    params:
        running_dir = Path(config['working_dir']) / 'trimming_ont' / 'summary'
    run:
        informs_filtering = SnakemakeUtils.load_object(Path(input.INFORMS_trimming))
        summary_data = [
            ('trim_ont_reads_in', informs_filtering['nb_seqs_in']),
            ('trim_ont_reads_out', informs_filtering['nb_seqs_out']),
            ('trim_ont_tool_version', informs_filtering['_name'])
        ]
        with open(output.TSV, 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')

rule trimming_ont_to_dict:
    """
    Combines the trimmed reads into a dictionary.
    """
    input:
        FASTQ = rules.trimming_ont_seqkit.output.FASTQ
    output:
        IO = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_DICT
    run:
        output_dict = {
            'SE': SnakemakeUtils.load_object(Path(input.FASTQ))
        }
        SnakemakeUtils.dump_object(output_dict,Path(output.IO))
