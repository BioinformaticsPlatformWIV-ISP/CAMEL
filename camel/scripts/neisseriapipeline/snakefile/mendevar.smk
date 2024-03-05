from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import sequence_typing
from camel.scripts.neisseriapipeline.snakefile import mendevar

rule mendevar_run:
    """
    Runs the MenDeVAR stand-alone tool.
    """
    input:
        TSV = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_TSV).format(
            scheme='bast', locus_type='peptide')
    output:
        TSV = Path(config['working_dir']) / 'mendevar' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'mendevar' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'mendevar',
        db = config.get('mendevar', {}).get('db')
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.tools.pipelines.neisseria.mendevar import MenDeVAR

        mendevar_ = MenDeVAR(Camel.get_instance())
        step = Step(str(rule), mendevar_, Camel.get_instance(), params.dir_)
        SnakemakeUtils.add_pickle_inputs(mendevar_,input)

        # Add database
        if params.db is None:
            raise RuntimeError('MenDeVAR database not specified in config')
        path_db = Path(params.db)
        if not path_db.exists():
            raise FileNotFoundError(f'MenDeVAR database not found: {path_db}')
        mendevar_.add_input_files({'DB': [ToolIOFile(path_db)]})

        # Run tool
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(mendevar_, output)

rule mendevar_report:
    """
    Creates an output report for the MenDeVAR analysis.
    """
    input:
        TSV = rules.mendevar_run.output.TSV,
        INFORMS_mendevar = rules.mendevar_run.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / mendevar.OUTPUT_MENDEVAR_REPORT
    params:
        dir_ = Path(config['working_dir']) / 'mendevar' / 'report'
    run:
        from camel.app.tools.pipelines.neisseria.mendevarreporter import MenDeVARReporter
        reporter = MenDeVARReporter(Camel.get_instance())
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule mendevar_report_empty:
    """
    Creates an empty MenDeVAR report when the analysis is disabled.
    """
    output:
        HTML = Path(config['working_dir']) / mendevar.OUTPUT_MENDEVAR_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('MenDeVAR', Path(output.HTML), 3)

rule mendevar_create_summary:
    """
    Creates the tabular summary output for the MenDeVAR assay.
    """
    input:
        INFORMS_mendevar = rules.mendevar_run.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / mendevar.OUTPUT_MENDEVAR_SUMMARY
    run:
        # Collect informs
        informs = SnakemakeUtils.load_object(Path(input.INFORMS_mendevar))

        # Create TSV output
        with open(output.TSV, 'w') as handle:
            handle.write(f"mendevar_bexsero_status\t{informs['bexsero_status']}")
            handle.write('\n')
            handle.write(f"mendevar_trumenba_status\t{informs['trumenba_status']}")
            handle.write('\n')
