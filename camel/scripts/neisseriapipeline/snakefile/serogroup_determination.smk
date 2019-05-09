from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.sequence_typing import OUTPUT_TYPING_HITS
from camel.scripts.neisseriapipeline.snakefile.serogroup_determination import OUTPUT_SEROGROUP_DETERMINATION_SUMMARY

camel = Camel.get_instance()


rule Serogroup_determination_analysis:
    """
    This rule is used to determine the serogroup based on the detected genes.
    """
    input:
        hits=expand(os.path.join(config['working_dir'], OUTPUT_TYPING_HITS),
                    locus_type='DNA',
                    detection_method=config['detection_method'],
                    scheme=sorted([k for k in config['sequence_typing'].keys() if k.startswith('serogroup')]))
    output:
        INFORMS=os.path.join(config['working_dir'], 'serogroup_determination', 'informs.io')
    params:
        serogroups=sorted([k for k in config['sequence_typing'].keys() if k.startswith('serogroup')]),
        working_dir=os.path.join(config['working_dir'], 'serogroup_determination')
    run:
        from camel.app.tools.pipelines.neisseria.serogroupdetermination import SerogroupDetermination
        serogroup_determination = SerogroupDetermination(camel)
        for hits_output, serogroup in zip(input.hits, params.serogroups):
            serogroup_determination.add_input_files({f'hits_{serogroup}': SnakemakeUtils.load_object(hits_output)})
        step = Step(rule, serogroup_determination, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(serogroup_determination, output)


rule Serogroup_determination_report:
    """
    This rule is used to determine the serogroup based on the detected genes.
    """
    input:
        INFORMS_analysis=os.path.join(config['working_dir'], 'serogroup_determination', 'informs.io')
    output:
        VAL_HTML=os.path.join(config['working_dir'], 'serogroup_determination', 'html.io')
    params:
        working_dir=os.path.join(config['working_dir'], 'serogroup_determination')
    run:
        from camel.app.tools.pipelines.neisseria.serogroupdeterminationreporter import SerogroupDeterminationReporter
        reporter = SerogroupDeterminationReporter(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)


rule Serogroup_determination_report_empty:
    """
    This rule is used to create an empty report for the serogroup determination
    """
    output:
        VAL_HTML=os.path.join(config['working_dir'], 'serogroup_determination', 'html-empty.io')
    params:
        working_dir=os.path.join(config['working_dir'], 'serogroup_determination')
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Serogroup determination', output.VAL_HTML)


rule Serogroup_determination_summary:
    """
    Collects the summary information for the serogroup determination.
    """
    input:
        INFORMS_analysis=os.path.join(config['working_dir'], 'serogroup_determination', 'informs.io')
    output:
        os.path.join(config['working_dir'], OUTPUT_SEROGROUP_DETERMINATION_SUMMARY)
    run:
        informs = SnakemakeUtils.load_object(input.INFORMS_analysis)
        with open(output[0], 'w') as handle:
            for k, v in [('detected_serogroup', informs['detected_serogroup']),
                         ('serogroup_nb_hits', informs['serogroups_sorted'][0]['nb_hits']),
                         ('serogroup_nb_hits_perfect', informs['serogroups_sorted'][0]['nb_hits_perfect']),
                         ('serogroup_total_loci', informs['serogroups_sorted'][0]['nb_loci_total'])]:
                handle.write('\t'.join([k, str(v)]))
                handle.write('\n')
