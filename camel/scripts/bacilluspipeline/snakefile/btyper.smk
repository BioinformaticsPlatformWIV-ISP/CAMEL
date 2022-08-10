from pathlib import Path

from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.bacilluspipeline.snakefile import btyper as bt

rule btyper_run:
    """
    Runs btyper on assembled contigs.
    """
    input:
        FASTA = Path(config['working_dir']) / bt.INPUT_BTYPER_FASTA
    output:
        TSV = Path(config['working_dir']) / bt.OUTPUT_VAL_BTYPER,
        INFORMS = Path(config['working_dir']) / bt.OUTPUT_INFORMS_BTYPER
    params:
        running_dir = Path(config['working_dir']) / 'btyper'
    run:
        from camel.app.tools.btyper.btyper import BTyper
        btyper = BTyper(camel)
        btyper.update_parameters(output_dir=params.running_dir)
        SnakemakeUtils.add_pickle_inputs(btyper,input)
        step = Step(rule,btyper,camel,params.running_dir,config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(btyper,output)

rule btyper_report:
    """
    Creates the report section for btyper.
    """
    input:
        TSV = rules.btyper_run.output.TSV,
        INFORMS_btyper = rules.btyper_run.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / bt.OUTPUT_BTYPER_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'btyper',
        sample_name= config['sample_name']
    run:
        from camel.app.tools.btyper.btyperreporter import BTyperReporter
        btyper_reporter = BTyperReporter(camel)
        SnakemakeUtils.add_pickle_inputs(btyper_reporter,input)
        btyper_reporter.update_parameters(sample_name=params.sample_name)
        step = Step(rule,btyper_reporter,camel,params.running_dir,config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(btyper_reporter,output)

rule btyper_report_empty:
    """
    Creates an empty HTML report for the BTyper analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / bt.OUTPUT_BTYPER_REPORT_EMPTY
    params:
        running_dir = Path(config['working_dir']) / 'btyper'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.btyper.btyperreporter import BTyperReporter
        SnakePipelineUtils.create_empty_report_section(BTyperReporter.TITLE,Path(output.VAL_HTML))

rule btyper_dump_summary_info:
    """
    Dumps the summary information for the BTyper workflow in tabular format.
    """
    input:
        INFORMS = Path(config['working_dir']) / rules.btyper_run.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / bt.OUTPUT_BTYPER_SUMMARY
    run:
        import json
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        data = []
        with open(output.TSV, 'w') as handle:
            handle.write('{}\t{}\n'.format('btyper_typing_results', json.dumps(data)))
