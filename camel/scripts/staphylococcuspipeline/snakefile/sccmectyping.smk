from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import gene_detection
from camel.scripts.staphylococcuspipeline.snakefile import sccmectyping


rule sccmec_typing_run:
    """
    Determines the SCCmec information based on the detected BLAST hits.
    """
    input:
        YAML = config['sccmec_typing']['profiles'],
        VAL_HITS = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_ALL_HITS).format(db='sccmec_genes')
    output:
        INFORMS = Path(config['working_dir']) / 'sccmec_typing' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'sccmec_typing'
    run:
        from camel.app.tools.sccmectyping.sccmectyping import SCCmecTyping
        from camel.app.io.tooliofile import ToolIOFile
        sccmec_typing = SCCmecTyping(Camel.get_instance())
        sccmec_typing.add_input_files({'YML': [ToolIOFile(input.YAML)]})
        SnakemakeUtils.add_pickle_input(sccmec_typing, 'VAL_HITS', input.VAL_HITS)
        step = Step(rule, sccmec_typing, Camel.get_instance(), params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sccmec_typing, output)

rule sccmec_typing_report:
    """
    Creates a report section for the SCCmec typing.
    """
    input:
        INFORMS = rules.sccmec_typing_run.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / sccmectyping.OUTPUT_SCCMEC_TYPING_REPORT
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.components.html.htmlreportsection import HtmlReportSection
        section = HtmlReportSection('SCC<i>mec</i> type', 3)
        informs = SnakemakeUtils.load_object(input.INFORMS)
        section.add_table(
            [[f'{key}:', value] for key, value in informs.items() if not key.startswith('_')],
            table_attributes=[('class', 'information')])
        SnakemakeUtils.dump_object([ToolIOValue(section)], output.HTML)

rule sccmec_typing_report_empty:
    """
    Creates an empty report when SCCmec typing is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / sccmectyping.OUTPUT_SCCMEC_TYPING_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('SCC<i>mec</i> type', output.VAL_HTML)

rule sccmec_typing_summary:
    """
    Creates the summary output for the SCCmec workflow.
    """
    input:
        INFORMS = rules.sccmec_typing_run.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / sccmectyping.OUTPUT_SCCMEC_TYPING_SUMMARY
    run:
        informs = SnakemakeUtils.load_object(input.INFORMS)
        with open(output.TSV, 'w') as handle:
            for complex_ in informs['complexes']:
                handle.write('\t'.join([complex_['key'], complex_['value']]))
                handle.write('\n')
