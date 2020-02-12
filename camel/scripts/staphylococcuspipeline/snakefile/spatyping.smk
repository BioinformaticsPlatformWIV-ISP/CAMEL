from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import assembly_spades
from camel.scripts.staphylococcuspipeline.snakefile import spatyping as spatyping_workflow

rule spatyping_blastn:
    """
    Runs BLASTN to align the sequences against the spa typing database.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA,
        DB_BLAST = config['spa_typing']['db']
    output:
        TSV = Path(config['working_dir']) / 'spatyping' / 'blastn' / 'blast_hits.tsv'
    params:
        running_dir = Path(config['working_dir']) / 'spatyping' / 'blastn'
    run:
        from camel.app.tools.blast.blastn import Blastn
        from camel.app.tools.spatyping.spatyping import SpaTyping
        from camel.app.io.tooliofile import ToolIOFile
        blastn = Blastn(Camel.get_instance())
        SnakemakeUtils.add_pickle_input(blastn, 'FASTA', input.FASTA)
        blastn.add_input_files({'DB_BLAST': [ToolIOFile(input.DB_BLAST)]})
        blastn.update_parameters(output_format=SpaTyping.BLASTN_OUTPUT_FORMAT)
        step = Step(rule, blastn, Camel.get_instance(), params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blastn, output)

rule spatyping_run:
    """
    Determines the spa-type based on the blastn output.
    """
    input:
        TSV = rules.spatyping_blastn.output.TSV,
        CSV_profiles = config['spa_typing']['profiles']
    output:
        VAL_hits = Path(config['working_dir']) / 'spatyping' / 'detection' / 'hits.io',
        INFORMS = Path(config['working_dir']) / 'spatyping' / 'detection' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'spatyping' / 'detection'
    run:
        from camel.app.tools.spatyping.spatyping import SpaTyping
        spatyping = SpaTyping(Camel.get_instance())
        SnakemakeUtils.add_pickle_input(spatyping, 'TSV', input.TSV)
        spatyping.add_input_files({'CSV_profiles': [ToolIOFile(input.CSV_profiles)]})
        step = Step(rule, spatyping, Camel.get_instance(), params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spatyping, output)

rule spatyping_report:
    """
    Creates a report for the spa-typing assay.
    """
    input:
        VAL_hits = rules.spatyping_run.output.VAL_hits,
        INFORMS_spa_typing = rules.spatyping_run.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / spatyping_workflow.OUTPUT_SPATYPING_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'spatyping' / 'report'
    run:
        from camel.app.tools.spatyping.spatypingreporter import SpaTypingReporter
        reporter = SpaTypingReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, Camel.get_instance(), params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule spatyping_report_empty:
    """
    Creates an empty report when spatyping is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / spatyping_workflow.OUTPUT_SPATYPING_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('<i>spa</i> typing', output.VAL_HTML)
