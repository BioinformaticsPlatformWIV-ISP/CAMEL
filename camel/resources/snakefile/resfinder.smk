from pathlib import Path

from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import resfinder as rf

rule resfinder_run:
    """
    Runs resfinder on assembled contigs.
    """
    input:
        FASTA = Path(config['working_dir']) / rf.INPUT_RESFINDER_FASTA,
    output:
        TSV_pheno_general = Path(config['working_dir']) / rf.OUTPUT_RESFINDER_PHENO,
        TSV_genes = Path(config['working_dir']) / rf.OUTPUT_RESFINDER_GENE,
        INFORMS = Path(config['working_dir']) / rf.OUTPUT_RESFINDER_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'resfinder',
        min_cov = config['resfinder']['min_cov'],
        threshold = config['resfinder']['threshold'],
        acq_overlap = config['resfinder']['acq_overlap'],
        db = '/db/resfinder4'
    run:
        from camel.app.tools.resfinder.resfinder import ResFinder
        resfinder = ResFinder(camel)
        resfinder.update_parameters(output_path=str(params.running_dir), acquired=True,
            min_cov=params.min_cov, threshold=params.threshold, acq_overlap=params.acq_overlap)
        SnakemakeUtils.add_pickle_inputs(resfinder,input)
        resfinder.add_input_files({'DIR': [ToolIODirectory(Path(params.db))]})
        step = Step(rule,resfinder,camel,params.running_dir,config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(resfinder,output)

rule resfinder_report:
    """
    Creates the report section for resfinder.
    """
    input:
        TSV_genes = rules.resfinder_run.output.TSV_genes,
        TSV_pheno_general = rules.resfinder_run.output.TSV_pheno_general,
        INFORMS_resfinder = rules.resfinder_run.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / rf.OUTPUT_RESFINDER_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'resfinder',
        sample_name= config['sample_name']
    run:
        from camel.app.tools.resfinder.resfinderreporter import ResFinderReporter
        resfinder_reporter = ResFinderReporter(camel)
        SnakemakeUtils.add_pickle_inputs(resfinder_reporter,input)
        resfinder_reporter.update_parameters(sample_name=params.sample_name)
        step = Step(rule,resfinder_reporter,camel,params.running_dir,config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(resfinder_reporter,output)

rule resfinder_report_empty:
    """
    Creates an empty HTML report for the ResFinder analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / rf.OUTPUT_RESFINDER_REPORT_EMPTY
    params:
        running_dir = Path(config['working_dir']) / 'resfinder'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.resfinder.resfinderreporter import ResFinderReporter
        SnakePipelineUtils.create_empty_report_section(ResFinderReporter.TITLE,Path(output.VAL_HTML))

rule resfinder_dump_summary_info:
    """
    Dumps the summary information for the ResFinder workflow in tabular format.
    """
    input:
        INFORMS = Path(config['working_dir']) / rules.resfinder_run.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / rf.OUTPUT_RESFINDER_SUMMARY
    run:
        import json
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        data = []
        with open(output.TSV, 'w') as handle:
            handle.write('{}\t{}\n'.format('resfinder_typing_results', json.dumps(data)))
