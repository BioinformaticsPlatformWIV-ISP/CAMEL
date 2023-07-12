from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.resources.snakefile import sequence_typing, assembly_spades
from camel.scripts.neisseriapipeline.snakefile import serogroup_determination


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
        INFORMS = Path(config['working_dir']) / 'serogroup_determination' / 'legacy' / 'informs.io'
    params:
        serogroups = sorted([k for k in config['sequence_typing'].keys() if k.startswith('serogroup')]),
        working_dir = Path(config['working_dir']) / 'serogroup_determination'
    run:
        from camel.app.tools.pipelines.neisseria.serogroupdetermination import SerogroupDetermination
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        detector = SerogroupDetermination(Camel.get_instance())
        for hits_output, serogroup in zip(input.hits, params.serogroups):
            detector.add_input_files({f'hits_{serogroup}': SnakemakeUtils.load_object(Path(hits_output))})
        step = Step(str(rule), detector, Camel.get_instance(), params.working_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(detector, output)

rule serogroup_determination_report:
    """
    This rule is used to determine the serogroup based on the detected genes.
    """
    input:
        INFORMS_analysis = rules.serogroup_determination_analysis.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_LEGACY_REPORT
    params:
        working_dir = Path(config['working_dir']) / 'serogroup_determination'
    run:
        from camel.app.tools.pipelines.neisseria.serogroupdeterminationreporter import SerogroupDeterminationReporter
        reporter = SerogroupDeterminationReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, Camel.get_instance(), params.working_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule serogroup_determination_report_empty:
    """
    This rule is used to create an empty report for the serogroup determination
    """
    output:
        VAL_HTML = Path(config['working_dir']) / serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_LEGACY_REPORT_EMPTY
    params:
        working_dir = Path(config['working_dir']) / 'serogroup_determination'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Serogroup determination (legacy)', Path(output.VAL_HTML))

rule serogroup_capsule_tool:
    """
    Runs the serogroup capsule tool.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA
    output:
        TSV = Path(config['working_dir']) / 'serogroup_determination' / 'capsule' / 'tsv.io',
        JSON = Path(config['working_dir']) / 'serogroup_determination' / 'capsule' / 'json.io',
        INFORMS = Path(config['working_dir']) / 'serogroup_determination' / 'capsule' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'serogroup_determination' / 'capsule'
    run:
        from camel.app.tools.pipelines.neisseria.characterizeneisseriacapsule import CharacterizeNeisseriaCapsule
        capsule_typer = CharacterizeNeisseriaCapsule(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(capsule_typer, input)
        step = Step(str(rule), capsule_typer, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(capsule_typer, output)

rule serogroup_capsule_tool_report:
    """
    Creates an output report for the capsule tool.
    """
    input:
        TSV = rules.serogroup_capsule_tool.output.TSV,
        JSON = rules.serogroup_capsule_tool.output.JSON,
        INFORMS_detector = rules.serogroup_capsule_tool.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_REPORT
    params:
        working_dir = Path(config['working_dir']) / 'serogroup_determination' / 'capsule'
    run:
        from camel.app.tools.pipelines.neisseria.characterizeneisseriacapsulereporter import \
            CharacterizeNeisseriaCapsuleReporter
        reporter = CharacterizeNeisseriaCapsuleReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, Camel.get_instance(), params.working_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule serogroup_capsule_tool_report_empty:
    """
    Creates an empty report when the analysis is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_REPORT_EMPTY
    params:
        working_dir = Path(config['working_dir']) / 'serogroup_determination'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Capsule characterization', Path(output.VAL_HTML))

rule serogroup_determination_summary:
    """
    Collects the summary information for the serogroup determination.
    """
    input:
        INFORMS_capsule = rules.serogroup_capsule_tool.output.INFORMS,
        INFORMS_legacy = rules.serogroup_determination_analysis.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_SUMMARY
    run:
        rows_out = []

        # Capsule tool
        informs_capsule = SnakemakeUtils.load_object(Path(input.INFORMS_capsule))
        rows_out.extend([
            ('serogroup_capsule', informs_capsule['detected_serogroup']),
            ('serogroup_capsule_genes', informs_capsule['genes_present'])
        ])

        # Legacy detection
        informs_legacy = SnakemakeUtils.load_object(Path(input.INFORMS_legacy))
        rows_out.extend([
            ('serogroup_legacy', informs_legacy['detected_serogroup']),
            ('serogroup_legacy_nb_hits', informs_legacy['serogroups_sorted'][0]['nb_hits']),
            ('serogroup_legacy_nb_hits_perfect', informs_legacy['serogroups_sorted'][0]['nb_hits_perfect']),
            ('serogroup_legacy_total_loci', informs_legacy['serogroups_sorted'][0]['nb_loci_total'])
        ])

        # Write output file
        with open(output.TSV, 'w') as handle:
            for k, v in rows_out:
                handle.write('\t'.join([k, str(v)]))
                handle.write('\n')
