from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import mykrobe

rule mykrobe_run:
    """
    Runs the Mykrobe assay.
    """
    input:
        IO = mykrobe.get_input(config)
    output:
        CSV = Path(config['working_dir']) / 'mykrobe' / 'csv.io',
        INFORMS = Path(config['working_dir']) / 'mykrobe' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'mykrobe',
        species = Path(config['mykrobe']['species']),
        db_dir = Path(config['mykrobe']['db']),
        input_type = config['input_type']
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.mykrobe.mykrobe import Mykrobe

        typer = Mykrobe(Camel.get_instance())

        # Extract FASTQ paths to add them as Mykrobe input
        if params.input_type == 'illumina':
            fq_in = FastqInput.from_fq_dict(Path(input.IO), params.input_type)
            typer.add_input_files({'FASTQ_PE': fq_in.pe})
        if params.input_type == 'fasta':
            SnakemakeUtils.add_pickle_input(typer, 'FASTA', Path(input.IO))
        typer.add_input_files({
            'DIR': [ToolIODirectory(params.db_dir)],
            'SPECIES': [ToolIOValue(params.species)]
        })

        # Run tool
        step = Step(str(rule), typer, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(typer, output)


rule mykrobe_report:
    """
    Creates an output report for the Mykrobe analysis.
    """
    input:
        CSV = rules.mykrobe_run.output.CSV,
        INFORMS_mykrobe = rules.mykrobe_run.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / mykrobe.OUTPUT_MYKROBE_REPORT
    params:
        dir_ = Path(config['working_dir']) / 'mykrobe' / 'report',
        skip_amr = config['mykrobe'].get('skip_amr', False),
        title = config['mykrobe'].get('title', 'Lineage information')
    run:
        from camel.app.tools.mykrobe.mykrobereporter import MykrobeReporter

        reporter = MykrobeReporter(Camel.get_instance())
        if params.skip_amr:
            reporter.add_input_files({
                'SKIP_AMR': [ToolIOValue(params.skip_amr)],
                'custom_header': [ToolIOValue(params.title)]
            })
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule mykrobe_report_empty:
    """
    Creates an empty Mykrobe report when the analysis is disabled.
    """
    output:
        HTML = Path(config['working_dir']) / mykrobe.OUTPUT_MYKROBE_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Mykrobe', Path(output.HTML), 3)

rule mykrobe_create_summary:
    """
    Creates the tabular summary output for the Mykrobe assay.
    """
    input:
        INFORMS_mykrobe = rules.mykrobe_run.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / mykrobe.OUTPUT_MYKROBE_SUMMARY
    params:
        skip_amr = config['mykrobe'].get('skip_amr',False)
    run:
        # Collect informs
        informs = SnakemakeUtils.load_object(Path(input.INFORMS_mykrobe))

        # Create TSV output
        with open(output.TSV, 'w') as handle:
            handle.write(f"mykrobe_phylo_group\t{informs['phylo_group']}")
            handle.write('\n')
            handle.write(f"mykrobe_species\t{informs['species']}")
            handle.write('\n')
            handle.write(f"mykrobe_lineage\t{informs['lineage']}")
            handle.write('\n')
            if str(params.skip_amr) == 'False':
                handle.write(f"mykrobe_drug_susceptibility\t{informs['drug_susceptibility']}")
                handle.write('\n')
            handle.write(f"mykrobe_tool_version\t{informs['_name']}")
            handle.write('\n')
            handle.write(f"mykrobe_db_version\t{informs['db_version']}")
            handle.write('\n')
