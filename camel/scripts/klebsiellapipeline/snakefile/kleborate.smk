from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import assembly_spades

rule kleborate_run:
    """
    Runs the Kleborate tool.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA
    output:
        TSV = Path(config['working_dir']) / 'kleborate' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'kleborate' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'kleborate'
    run:
        from camel.app.tools.kleborate.kleborate import Kleborate
        kleborate = Kleborate(Camel.get_instance())
        kleborate.update_parameters(all=True)
        SnakemakeUtils.add_pickle_inputs(kleborate, input)
        step = Step(str(rule), kleborate, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(kleborate, output)

rule kleborate_reporter:
    """
    Creates the output report for Kleborate.
    """
    input:
        TSV = rules.kleborate_run.output.TSV,
        INFORMS_kleborate = rules.kleborate_run.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / 'kleborate' / 'html.io'
    params:
        dir_ = Path(config['working_dir']) / 'kleborate'
    run:
        from camel.app.tools.kleborate.kleboratereporter import KleborateReporter
        reporter = KleborateReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule kleborate_report_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'kleborate' / 'html-empty.io'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Kleborate', Path(output.VAL_HTML))

rule kleborate_create_summary:
    """
    Creates a tabular summary output for Kleborate.
    """
    input:
        INFORMS = rules.kleborate_run.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / 'kleborate' / 'summary_kleborate.tsv'
    params:
        keys_kept = ['ST', 'Yersiniabactin', 'YbST', 'wzi', 'K_locus', 'virulence_score']
    run:
        informs_in = SnakemakeUtils.load_object(Path(input.INFORMS))
        with open(output.TSV, 'w') as handle:
            for key in params.keys_kept:
                handle.write('\t'.join([
                    f'kleborate_{key}', informs_in[key]
                ]))
                handle.write('\n')
