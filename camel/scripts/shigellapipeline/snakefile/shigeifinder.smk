from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.resources.snakefile import assembly
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.shigellapipeline.snakefile import shigeifinder

rule shigeifinder_run:
    """
    Runs the ShigEiFinder assay.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA
    output:
        TSV = Path(config['working_dir']) / 'shigeifinder' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'shigeifinder' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'shigeifinder'
    run:
        from camel.app.tools.pipelines.shigella.shigeifinder import ShigEiFinder

        shigeifinder_ = ShigEiFinder(Camel.get_instance())
        step = Step(str(rule), shigeifinder_, Camel.get_instance(), params.dir_)
        SnakemakeUtils.add_pickle_inputs(shigeifinder_,input)

        # Run tool
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(shigeifinder_, output)

rule shigeifinder_report:
    """
    Creates an output report for the ShigEiFinder analysis.
    """
    input:
        TSV = rules.shigeifinder_run.output.TSV,
        INFORMS_shigeifinder = rules.shigeifinder_run.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / shigeifinder.OUTPUT_SHIGEIFINDER_REPORT
    params:
        dir_ = Path(config['working_dir']) / 'shigeifinder' / 'report'
    run:
        from camel.app.tools.pipelines.shigella.shigeifinderreporter import ShigEiFinderReporter
        reporter = ShigEiFinderReporter(Camel.get_instance())
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule shigeifinder_report_empty:
    """
    Creates an empty ShigEiFinder report when the analysis is disabled.
    """
    output:
        HTML = Path(config['working_dir']) / shigeifinder.OUTPUT_SHIGEIFINDER_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('ShigEiFinder', Path(output.HTML), 3)

rule shigeifinder_create_summary:
    """
    Creates the tabular summary output for the ShigEiFinder assay.
    """
    input:
        INFORMS_shigeifinder = rules.shigeifinder_run.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / shigeifinder.OUTPUT_SHIGEIFINDER_SUMMARY
    run:
        # Collect informs
        informs = SnakemakeUtils.load_object(Path(input.INFORMS_shigeifinder))

        # Create TSV output
        with open(output.TSV, 'w') as handle:
            handle.write(f"shigeifinder_serotype\t{informs['species']}")
            handle.write('\n')
