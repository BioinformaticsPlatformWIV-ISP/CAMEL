from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.shigellapipeline.snakefile import shigatyper

rule shigatyper_run:
    """
    Runs the ShigaTyper assay.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        TSV = Path(config['working_dir']) / 'shigatyper' / 'tsv.io',
        TSV_HITS = Path(config['working_dir']) / 'shigatyper' / 'tsv_hits.io',
        INFORMS = Path(config['working_dir']) / 'shigatyper' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'shigatyper',
        input_type = config['input_type']
    run:
        from camel.app.tools.pipelines.shigella.shigatyper import ShigaTyper
        from camel.app.components.workflows.utils.fastqinput import FastqInput

        typer = ShigaTyper(Camel.get_instance())

        # Extract FASTQ paths to add them as ShigaTyper input
        fq_in = FastqInput.from_fq_dict(Path(input.IO),'illumina')
        typer.add_input_files({'FASTQ_PE': fq_in.pe})

        # Run tool
        step = Step(str(rule), typer, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(typer, output)

rule shigatyper_report:
    """
    Creates an output report for the ShigaTyper analysis.
    """
    input:
        TSV = rules.shigatyper_run.output.TSV,
        TSV_HITS = rules.shigatyper_run.output.TSV_HITS,
        INFORMS_shigatyper = rules.shigatyper_run.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / shigatyper.OUTPUT_SHIGATYPER_REPORT
    params:
        dir_ = Path(config['working_dir']) / 'shigatyper' / 'report'
    run:
        from camel.app.tools.pipelines.shigella.shigatyperreporter import ShigaTyperReporter
        reporter = ShigaTyperReporter(Camel.get_instance())
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule shigatyper_report_empty:
    """
    Creates an empty ShigaTyper report when the analysis is disabled.
    """
    output:
        HTML = Path(config['working_dir']) / shigatyper.OUTPUT_SHIGATYPER_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('ShigaTyper', Path(output.HTML), 3)

rule shigatyper_create_summary:
    """
    Creates the tabular summary output for the ShigaTyper assay.
    """
    input:
        INFORMS_shigatyper = rules.shigatyper_run.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / shigatyper.OUTPUT_SHIGATYPER_SUMMARY
    run:
        # Collect informs
        informs = SnakemakeUtils.load_object(Path(input.INFORMS_shigatyper))

        # Create TSV output
        with open(output.TSV, 'w') as handle:
            handle.write(f"shigatyper_prediction\t{informs['species']}")
            handle.write('\n')
            handle.write(f"shigatyper_tool_version\t{informs['_name']}")
            handle.write('\n')
