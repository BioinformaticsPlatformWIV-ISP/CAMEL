from pathlib import Path

from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import pointfinder, assembly_spades


rule pointfinder_run:
    """
    This rule executes PointFinder on the input sample.
    """
    input:
        FASTA = Path(config['working_dir'])/ assembly_spades.OUTPUT_ASSEMBLY_FASTA
    output:
        TSV = Path(config['working_dir']) / 'pointfinder' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / pointfinder.OUTPUT_POINTFINDER_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'pointfinder',
        db = config.get('pointfinder', {}).get('db', 'escherichia_coli')
    run:
        from camel.app.tools.pointfinder.pointfinder import PointFinder
        pointfinder_ = PointFinder(camel)
        SnakemakeUtils.add_pickle_inputs(pointfinder_, input)
        step = Step(rule, pointfinder_, camel, params.running_dir, config)
        pointfinder_.update_parameters(database=params.db)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(pointfinder_, output)

rule pointfinder_report:
    """
    This rule creates the report for the pointfinder assay.
    """
    input:
        TSV = rules.pointfinder_run.output.TSV,
        INFORMS_pointfinder = rules.pointfinder_run.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / pointfinder.OUTPUT_POINTFINDER_REPORT,
        INFORMS = Path(config['working_dir']) / 'pointfinder' / 'informs-report.io'
    params:
        running_dir = Path(config['working_dir']) / 'pointfinder',
        sample_name = config['sample_name']
    run:
        from camel.app.tools.pointfinder.pointfinderreporter import PointFinderReporter
        pointfinder_reporter = PointFinderReporter(camel)
        SnakemakeUtils.add_pickle_inputs(pointfinder_reporter, input)
        pointfinder_reporter.update_parameters(sample_name=params.sample_name)
        step = Step(rule, pointfinder_reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(pointfinder_reporter, output)

rule pointfinder_report_empty:
    """
    Creates an empty HTML report for the PointFinder analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / pointfinder.OUTPUT_POINTFINDER_REPORT_EMPTY
    params:
        running_dir = Path(config['working_dir']) / 'pointfinder'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.pointfinder.pointfinderreporter import PointFinderReporter
        SnakePipelineUtils.create_empty_report_section(PointFinderReporter.TITLE, Path(output.VAL_HTML))

rule pointfinder_dump_summary_info:
    """
    Dumps the summary information for the PointFinder workflow in tabular format.
    """
    input:
        INFORMS = Path(config['working_dir']) / 'pointfinder' / 'informs-report.io'
    output:
        TSV = Path(config['working_dir']) / pointfinder.OUTPUT_POINTFINDER_SUMMARY
    run:
        import json
        from camel.app.components.html.htmltablecell import HtmlTableCell
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        mutations = []
        for row in informs['mutations']:
            mutations.append([e if not isinstance(e, HtmlTableCell) else e.text for e in row])
        with open(output.TSV, 'w') as handle:
            handle.write('{}\t{}'.format('pointfinder_mutations', json.dumps(mutations)))
            handle.write('\n')
