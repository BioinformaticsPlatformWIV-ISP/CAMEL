from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.scripts.yersiniapipeline.snakefile import species_determination

rule species_determination_analysis:
    """
    This rule is used to determine the species based on the cgST.
    """
    input:
        TSV_profile_matches = Path(config['working_dir']) / 'species_determination' / 'input' / 'tsv_profile_matches.io',
        TSV_taxonomic = Path(config['working_dir']) / 'species_determination' / 'input' / 'tsv_taxonomic.io'
    output:
        INFORMS = Path(config['working_dir']) / species_determination.OUTPUT_SPECIES_DETERMINATION_INFORMS,
        TSV = Path(config['working_dir']) / species_determination.OUTPUT_SPECIES_DETERMINATION_TSV
    params:
        working_dir = Path(config['working_dir']) / 'species_determination'
    run:
        from camel.app.tools.pipelines.yersinia.speciesdetermination import SpeciesDetermination
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        detector = SpeciesDetermination(Camel.get_instance())
        detector.add_input_files({'profile_matches': SnakemakeUtils.load_object(Path(str(input.TSV_profile_matches))),
                                  'taxonomic_file': SnakemakeUtils.load_object(Path(str(input.TSV_taxonomic)))})
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

rule species_determination_create_summary:
    """
    Creates the tabular summary output for the species determination.
    """
    input:
        INFORMS = Path(config['working_dir']) / species_determination.OUTPUT_SPECIES_DETERMINATION_INFORMS
    output:
        TSV = Path(config['working_dir']) / species_determination.OUTPUT_SPECIES_DETERMINATION_SUMMARY
    run:
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))['best_match']
        with open(output.TSV, 'w') as handle:
            for k, v in informs.items():
                handle.write('\t'.join(["species_determination_"+k, str(v)]))
                handle.write('\n')