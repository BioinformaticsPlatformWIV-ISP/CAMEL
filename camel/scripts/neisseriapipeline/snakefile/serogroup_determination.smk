from pathlib import Path

from camel.app.camel import Camel
from camel.resources.snakefile import sequence_typing
from camel.scripts.neisseriapipeline.snakefile import serogroup_determination

camel = Camel.get_instance()


rule serogroup_determination_analysis:
    """
    This rule is used to determine the serogroup based on the detected genes.
    """
    input:
        hits = expand(str(Path(config['working_dir']) / sequence_typing.OUTPUT_TYPING_HITS),
         locus_type='DNA',
         detection_method = config['detection_method'],
         scheme = sorted([k for k in config['sequence_typing'].keys() if k.startswith('serogroup')]))
    output:
        INFORMS = Path(config['working_dir']) / 'serogroup_determination' / 'informs.io'
    params:
        serogroups = sorted([k for k in config['sequence_typing'].keys() if k.startswith('serogroup')]),
        working_dir = Path(config['working_dir']) / 'serogroup_determination'
    run:
        from camel.app.tools.pipelines.neisseria.serogroupdetermination import SerogroupDetermination
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        detector = SerogroupDetermination(camel)
        for hits_output, serogroup in zip(input.hits, params.serogroups):
            detector.add_input_files({f'hits_{serogroup}': SnakemakeUtils.load_object(Path(hits_output))})
        step = Step(rule, detector, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(detector, output)


rule serogroup_determination_report:
    """
    This rule is used to determine the serogroup based on the detected genes.
    """
    input:
        INFORMS_analysis = rules.serogroup_determination_analysis.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_REPORT
    params:
        working_dir = Path(config['working_dir']) / 'serogroup_determination'
    run:
        from camel.app.tools.pipelines.neisseria.serogroupdeterminationreporter import SerogroupDeterminationReporter
        reporter = SerogroupDeterminationReporter(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)


rule serogroup_determination_report_empty:
    """
    This rule is used to create an empty report for the serogroup determination
    """
    output:
        VAL_HTML = Path(config['working_dir']) / serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_REPORT_EMPTY
    params:
        working_dir = Path(config['working_dir']) / 'serogroup_determination'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Serogroup determination', Path(output.VAL_HTML))


rule serogroup_determination_summary:
    """
    Collects the summary information for the serogroup determination.
    """
    input:
        INFORMS_analysis = rules.serogroup_determination_analysis.output.INFORMS
    output:
        TSV =Path(config['working_dir']) / serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_SUMMARY
    run:
        informs = SnakemakeUtils.load_object(Path(input.INFORMS_analysis))
        with open(output.TSV, 'w') as handle:
            for k, v in [('detected_serogroup', informs['detected_serogroup']),
                         ('serogroup_nb_hits', informs['serogroups_sorted'][0]['nb_hits']),
                         ('serogroup_nb_hits_perfect', informs['serogroups_sorted'][0]['nb_hits_perfect']),
                         ('serogroup_total_loci', informs['serogroups_sorted'][0]['nb_loci_total'])]:
                handle.write('\t'.join([k, str(v)]))
                handle.write('\n')
