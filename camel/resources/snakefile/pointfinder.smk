import json
import os

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.assembly_spades import OUTPUT_ASSEMBLY_FASTA
from camel.resources.snakefile.pointfinder import OUTPUT_POINTFINDER_REPORT, OUTPUT_POINTFINDER_REPORT_EMPTY, \
    OUTPUT_POINTFINDER_SUMMARY

rule PointFinder_run:
    """
    This rule executes PointFinder on the input sample.
    """
    input:
        FASTA=os.path.join(config['working_dir'], OUTPUT_ASSEMBLY_FASTA)
    output:
        TSV=os.path.join(config['working_dir'], 'pointfinder', 'tsv.io'),
        INFORMS=os.path.join(config['working_dir'], 'pointfinder', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'pointfinder')
    run:
        from camel.app.tools.pointfinder.pointfinder import PointFinder
        pointfinder = PointFinder(camel)
        SnakemakeUtils.add_pickle_inputs(pointfinder, input)
        step = Step(rule, pointfinder, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(pointfinder, output)

rule PointFinder_report:
    """
    This rule creates the report for the pointfinder assay.
    """
    input:
        TSV=os.path.join(config['working_dir'], 'pointfinder', 'tsv.io'),
        INFORMS_pointfinder=os.path.join(config['working_dir'], 'pointfinder', 'informs.io')
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_POINTFINDER_REPORT),
        INFORMS=os.path.join(config['working_dir'], 'pointfinder', 'informs-report.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'pointfinder'),
        sample_name=config['sample_name']
    run:
        from camel.app.tools.pointfinder.pointfinderreporter import PointFinderReporter
        pointfinder_reporter = PointFinderReporter(camel)
        SnakemakeUtils.add_pickle_inputs(pointfinder_reporter, input)
        pointfinder_reporter.update_parameters(sample_name=params.sample_name)
        step = Step(rule, pointfinder_reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(pointfinder_reporter, output)

rule PointFinder_report_empty:
    """
    Creates an empty HTML report for the PointFinder analysis.
    """
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_POINTFINDER_REPORT_EMPTY)
    params:
        running_dir=os.path.join(config['working_dir'], 'pointfinder')
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.pointfinder.pointfinderreporter import PointFinderReporter
        SnakePipelineUtils.create_empty_report_section(PointFinderReporter.TITLE, output.VAL_HTML)

rule PointFinder_dump_summary_info:
    """
    Dumps the summary information for the PointFinder workflow in tabular format.
    """
    input:
        INFORMS=os.path.join(config['working_dir'], 'pointfinder', 'informs-report.io')
    output:
        os.path.join(config['working_dir'], OUTPUT_POINTFINDER_SUMMARY)
    run:
        informs = SnakemakeUtils.load_object(input.INFORMS)
        with open(output[0], 'w') as handle:
            handle.write('{}\t{}'.format('pointfinder_mutations', json.dumps(informs['mutations'])))
            handle.write('\n')
