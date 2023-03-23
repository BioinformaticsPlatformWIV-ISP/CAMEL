from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import assembly_spades


rule mobsuite_mob_recon:
    """
    Runs the MOB-recon tool.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA,
        DB = config['mob_suite']['db']
    output:
        TSV = Path(config['working_dir']) / 'mob_suite' / 'tsv.io',
        TSV_contigs = Path(config['working_dir']) / 'mob_suite' / 'tsv-contigs.io',
        FASTA = Path(config['working_dir']) / 'mob_suite' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / 'mob_suite' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'mob_suite'
    threads: 4
    run:
        from camel.app.tools.mobsuite.mobrecon import MOBRecon
        mob_recon = MOBRecon(Camel.get_instance())
        SnakemakeUtils.add_pickle_input(mob_recon, 'FASTA', Path(input.FASTA))
        mob_recon.add_input_files({'DB': [ToolIODirectory(Path(input.DB))]})
        mob_recon.update_parameters(num_threads=threads)
        step = Step(str(rule), mob_recon, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(mob_recon, output)

rule mobsuite_mob_recon_reporter:
    """
    Creates the output report for MOB-recon.
    """
    input:
        TSV = rules.mobsuite_mob_recon.output.TSV,
        TSV_contigs = rules.mobsuite_mob_recon.output.TSV_contigs,
        FASTA = rules.mobsuite_mob_recon.output.FASTA,
        INFORMS_mob_recon = rules.mobsuite_mob_recon.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / 'mob_suite' / 'html.io'
    params:
        dir_ = Path(config['working_dir']) / 'mob_suite'
    run:
        from camel.app.tools.mobsuite.mobreconreporter import MOBReconReporter
        reporter = MOBReconReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule mobsuite_mob_recon_report_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'mob_suite' / 'html-empty.io'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('MOB-recon', Path(output.VAL_HTML))
