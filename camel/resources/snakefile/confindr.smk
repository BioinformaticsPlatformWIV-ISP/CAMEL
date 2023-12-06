from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import confindr

rule confindr_run:
    """
    Runs ConFindr on the input FASTQ datasets.
    """
    input:
        IO_FASTQ = Path(config['working_dir']) / 'fq_dict.io'
    output:
        CSV = Path(config['working_dir']) / 'confindr' / 'csv.io',
        INFORMS = Path(config['working_dir']) / 'confindr' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'confindr',
        db = config.get('confindr', {}).get('db'),
        input_type = config['input_type']
    threads: 4
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.confindr.confindr import ConFindr

        confindr_ = ConFindr(Camel.get_instance())
        confindr_.update_parameters(rmlst=True, threads=threads)

        # Add database
        if params.db is None:
            raise RuntimeError('ConFindr database not specified in config')
        path_db = Path(params.db)
        if not path_db.exists():
            raise FileNotFoundError(f'ConFindr database not found: {path_db}')
        confindr_.update_parameters(rmlst=True, databases=str(path_db))

        # Add input files
        if params.input_type not in ('hybrid', 'illumina'):
            raise ValueError('ConFindr currently only support Illumina input')
        fq_in = FastqInput.from_fq_dict(Path(input.IO_FASTQ), 'illumina')
        confindr_.add_input_files({'FASTQ_PE': fq_in.pe})

        # Run the tool
        step = Step(str(rule), confindr_, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(confindr_, output)

rule confindr_report:
    """
    Creates an output report for the ConFindr analysis.
    """
    input:
        INFORMS_confindr = rules.confindr_run.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / confindr.OUTPUT_CONFINDR_REPORT
    params:
        dir_ = Path(config['working_dir']) / 'confindr' / 'report'
    run:
        from camel.app.tools.confindr.confindrreporter import ConFindrReporter
        reporter = ConFindrReporter(Camel.get_instance())
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule confindr_report_empty:
    """
    Creates an empty ConFindr report when the analysis is disabled.
    """
    output:
        HTML = Path(config['working_dir']) / confindr.OUTPUT_CONFINDR_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('ConFindr', Path(output.HTML), 2)

rule confindr_create_summary:
    """
    Creates the tabular summary output for the ConFindr assay.
    """
    input:
        INFORMS_confindr = rules.confindr_run.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / confindr.OUTPUT_CONFINDR_SUMMARY
    run:
        # Collect informs
        informs = SnakemakeUtils.load_object(Path(input.INFORMS_confindr))
        rows_out = [
            ('genus', informs['Genus']),
            ('nb_contam_snvs', informs['NumContamSNVs']),
            ('contam_status', informs['ContamStatus'])
        ]

        # Create TSV output
        with open(output.TSV, 'w') as handle:
            for key, value in rows_out:
                handle.write(f'confindr_{key}\t{value}')
                handle.write('\n')
