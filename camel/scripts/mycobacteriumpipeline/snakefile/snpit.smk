from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import variant_calling
from camel.scripts.mycobacteriumpipeline.snakefile import snpit as snpit_workflow


rule snpit_report:
    """
    Creates the report for the snpit assay.
    """
    input:
        VCF = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_UNFILTERED_VCF
    output:
        VAL_HTML = Path(config['working_dir']) / snpit_workflow.OUTPUT_SNPIT_REPORT,
        INFORMS = Path(config['working_dir']) / snpit_workflow.OUTPUT_SNPIT_INFORMS
    params:
        dir_ = Path(config['working_dir']) / 'snpit'
    run:
        from camel.app.tools.snpit.snpit import Snpit
        snpit = Snpit(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(snpit, input)
        step = Step(str(rule), snpit, Camel.get_instance(), Path(params.dir_))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(snpit, output)

rule snpit_report_empty:
    """
    Creates an empty report when the snpit assay is disabled.
    """
    output:
        HTML = Path(config['working_dir']) / snpit_workflow.OUTPUT_SNPIT_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Snpit', Path(output.HTML), 2)

rule snpit_export_summary:
    """
    Exports the summary information for the snpit assay.
    """
    input:
        INFORMS = rules.snpit_report.output.INFORMS
    output:
        TSV = Path(config['working_dir'], snpit_workflow.OUTPUT_SNPIT_SUMMARY)
    run:
        snpit_informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        with open(output.TSV, 'w') as handle:
            for key in ('species', 'lineage', 'sublineage', 'percent_matched'):
                handle.write('\t'.join([f'snpit_{key}', str(snpit_informs[key])]))
                handle.write('\n')
