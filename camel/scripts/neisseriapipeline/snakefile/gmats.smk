from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import sequence_typing
from camel.scripts.neisseriapipeline.snakefile import gmats

rule gmats_run:
    """
    Runs the gMATS assay.
    """
    input:
        TSV = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_TSV).format(
            scheme='bast', locus_type='peptide')
    output:
        TSV = Path(config['working_dir']) / 'gmats' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'gmats' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'gmats',
        db = config.get('gmats', {}).get('db')
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.tools.pipelines.neisseria.gmats import GMats

        gmats_ = GMats(Camel.get_instance())
        step = Step(str(rule), gmats_, Camel.get_instance(), params.dir_)
        SnakemakeUtils.add_pickle_inputs(gmats_,input)

        # Add database
        if params.db is None:
            raise RuntimeError('gMATS database not specified in config')
        path_db = Path(params.db)
        if not path_db.exists():
            raise FileNotFoundError(f'gMATS database not found: {path_db}')
        gmats_.add_input_files({'DB': [ToolIOFile(path_db)]})

        # Run tool
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(gmats_, output)

rule gmats_report:
    """
    Creates an output report for the gMATS analysis.
    """
    input:
        TSV = rules.gmats_run.output.TSV,
        INFORMS_gmats = rules.gmats_run.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / gmats.OUTPUT_GMATS_REPORT
    params:
        dir_ = Path(config['working_dir']) / 'gmats' / 'report'
    run:
        from camel.app.tools.pipelines.neisseria.gmatsreporter import GMatsReporter
        reporter = GMatsReporter(Camel.get_instance())
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule gmats_report_empty:
    """
    Creates an empty gMATS report when the analysis is disabled.
    """
    output:
        HTML = Path(config['working_dir']) / gmats.OUTPUT_GMATS_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('gMATS', Path(output.HTML), 3)

rule gmats_create_summary:
    """
    Creates the tabular summary output for the gMATS assay.
    """
    input:
        INFORMS_gmats = rules.gmats_run.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / gmats.OUTPUT_GMATS_SUMMARY
    run:
        # Collect informs
        informs = SnakemakeUtils.load_object(Path(input.INFORMS_gmats))

        # Create TSV output
        with open(output.TSV, 'w') as handle:
            handle.write(f"gmats_status\t{informs['gMATS_status']}")
            handle.write('\n')
