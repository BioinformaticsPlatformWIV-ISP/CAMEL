from pathlib import Path
import shutil

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import trimming_ont

camel = Camel.get_instance()

rule trimming_ont_pickle_fastq_input:
    """
    Creates a pickle for the fastq input files.
    """
    output:
        FASTQ_SE = Path(config['working_dir']) / 'trimming_ont' / 'input' / 'fastq-se.io'
    params:
        config_input = config['input']
    run:
        from camel.app.io.tooliofile import ToolIOFile
        if trimming_ont.INPUT_ONT_FASTQ.exists():
            shutil.copyfile(trimming_ont.INPUT_ONT_FASTQ, Path(output.FASTQ_SE))
        else:
            fastq_se_in = Path([params.config_input[fastq_key][0]['path'] for fastq_key in ['fastq_se', 'fastq'] if fastq_key in params.config_input][0])
            SnakemakeUtils.dump_object([ToolIOFile(fastq_se_in)], Path(output.FASTQ_SE))


rule trimming_ont_nanoplot_pre:
    """
    Creates NanoPlot reports for the raw reads. 
    """
    input:
        FASTQ = rules.trimming_ont_pickle_fastq_input.output.FASTQ_SE
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


rule trimming_ont_filtlong:
    """
    Read trimming using filtlong.
    """
    input:
        FASTQ = rules.trimming_ont_pickle_fastq_input.output.FASTQ_SE
    output:
        FASTQ = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_READS,
        INFORMS = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_INFORMS
    threads: 4
    priority: 1
    params:
        running_dir = Path(config['working_dir']) / 'trimming_ont' / 'filtlong',
        sample_name = config.get('sample_name', 'reads')
    run:
        from camel.app.tools.filtlong.filtlong import Filtlong
        filtlong = Filtlong(camel)
        SnakemakeUtils.add_pickle_inputs(filtlong, input)
        step = Step(str(rule), filtlong, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(filtlong, output)


rule trimming_ont_nanoplot_post:
    """
    Creates NanoPlot reports of the trimmed reads.
    """
    input:
        FASTQ = rules.trimming_ont_filtlong.output.FASTQ
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
        FASTQ_PE = rules.trimming_ont_filtlong.output.FASTQ,
        INFORMS_trimming = rules.trimming_ont_filtlong.output.INFORMS,
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
        INFORMS_trimming = rules.trimming_ont_filtlong.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_SUMMARY
    params:
        running_dir = Path(config['working_dir']) / 'trimming_ont' / 'summary'
    run:
        filtlong_informs = SnakemakeUtils.load_object(Path(input.INFORMS_trimming))
        summary_data = [
            ('nb_reads_in', filtlong_informs['nb_reads_in']),
            ('nb_reads_out', filtlong_informs['nb_reads_out'])
        ]
        with open(output[0], 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')


rule trimming_ont_to_dict:
    """
    Combines the trimmed reads into a dictionary.
    """
    input:
        FASTQ = rules.trimming_ont_filtlong.output.FASTQ
    output:
        IO = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_DICT
    run:
        output_dict = {
            'SE': SnakemakeUtils.load_object(Path(input.FASTQ))
        }
        SnakemakeUtils.dump_object(output_dict,Path(output.IO))
