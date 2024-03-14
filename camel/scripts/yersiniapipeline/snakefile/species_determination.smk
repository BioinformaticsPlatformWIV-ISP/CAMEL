from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.scripts.yersiniapipeline.snakefile import species_determination
from camel.resources.snakefile import sequence_typing

rule species_determination_analysis:
    """
    This rule is used to determine the species based on the cgST.
    """
    input:
        TSV_ST = Path(config['working_dir']) / 'species_determination' / 'input' / 'tsv.io'
    output:
        INFORMS= Path(config['working_dir']) / species_determination.OUTPUT_SPECIES_DETERMINATION_INFORMS,
        TSV = Path(config['working_dir']) / species_determination.OUTPUT_SPECIES_DETERMINATION_TSV
    params:
        working_dir = Path(config['working_dir']) / 'species_determination'
    run:
        from camel.app.tools.pipelines.yersinia.speciesdetermination import SpeciesDetermination
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        detector = SpeciesDetermination(Camel.get_instance())
        detector.add_input_files({'profile_matches': SnakemakeUtils.load_object(Path(str(input.TSV_ST)))})
        step = Step(str(rule), detector, Camel.get_instance(), params.working_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(detector, output)

rule species_determination_report:
    """
    This rule is used to determine the species determination HTML report.
    """
    input:
        INFORMS_analysis = rules.species_determination_analysis.output.INFORMS,
        TSV_analysis = rules.species_determination_analysis.output.TSV
    output:
        VAL_HTML = Path(config['working_dir']) / species_determination.OUTPUT_SPECIES_DETERMINATION_REPORT
    params:
        working_dir = Path(config['working_dir']) / 'species_determination'
    run:
        from camel.app.tools.pipelines.yersinia.speciesdeterminationreporter import SpeciesDeterminationReporter
        reporter = SpeciesDeterminationReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, Camel.get_instance(), params.working_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule species_determination_report_empty:
    """
    This rule is used to create an empty report for the species determination.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / species_determination.OUTPUT_SPECIES_DETERMINATION_REPORT_EMPTY
    params:
        working_dir = Path(config['working_dir']) / 'species_determination'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Species determination', Path(output.VAL_HTML))
