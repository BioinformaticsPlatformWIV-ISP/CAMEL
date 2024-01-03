from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import assembly

rule bacmet_pickle_db:
    """
    Creates a pickle object for the database.
    """
    input:
        DB = config['bacmet']['db']
    output:
        DB = Path(config['working_dir']) / 'bacmet' / 'db.io'
    run:
        from camel.app.io.tooliodirectory import ToolIODirectory
        SnakemakeUtils.dump_object([ToolIODirectory(Path(input.DB))], Path(output.DB))

rule bacmet_prodigal:
    """
    Runs Prodigal to predict CDS.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA
    output:
        FASTA = Path(config['working_dir']) / 'bacmet' / 'prodigal' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / 'bacmet' / 'prodigal' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'bacmet' / 'prodigal'
    run:
        from camel.app.tools.prodigal.prodigal import Prodigal
        prodigal = Prodigal(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(prodigal, input)
        step = Step(str(rule), prodigal, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(prodigal, output)

rule bacmet_prodigal_report:
    """
    Creates an output report for Prodigal.
    """
    input:
        FASTA = rules.bacmet_prodigal.output.FASTA,
        INFORMS_prodigal = rules.bacmet_prodigal.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / 'bacmet' / 'prodigal' / 'report' / 'html.io'
    params:
        dir_ = Path(config['working_dir']) / 'bacmet' / 'prodigal' / 'report'
    run:
        from camel.app.tools.prodigal.prodigalreporter import ProdigalReporter
        reporter = ProdigalReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule bacmet_prodigal_report_empty:
    """
    Creates an empty output report for Prodigal.
    """
    output:
        HTML = Path(config['working_dir']) / 'bacmet' / 'prodigal' / 'report' / 'html-empty.io'
    params:
        dir_ = Path(config['working_dir']) / 'bacmet' / 'prodigal' / 'report'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Prodigal', Path(output.HTML))

rule bacmet_blastp:
    """
    Runs blastp to identify protein matches.
    """
    input:
        FASTA = rules.bacmet_prodigal.output.FASTA,
        DB = rules.bacmet_pickle_db.output.DB
    output:
        TSV = Path(config['working_dir']) / 'bacmet' / 'blastp' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'bacmet' / 'blastp' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'bacmet' / 'blastp',
        fmt = '6 pident sseqid sseq slen qseqid qstart qend'
    threads: 4
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.tools.blast.blastp import Blastp

        # Create & run tool
        blastp = Blastp(Camel.get_instance())
        blastp.update_parameters(output_format=f'"{params.fmt}"', threads=threads)
        SnakemakeUtils.add_pickle_input(blastp, 'FASTA', Path(input.FASTA))
        path_db = next(SnakemakeUtils.load_object(Path(input.DB))[0].path.glob('*.fasta'))
        blastp.add_input_files({'DB_BLAST': [ToolIOFile(path_db)]})
        step = Step(str(rule), blastp, Camel.get_instance(), params.dir_)
        step.run_step()

        # Dump output
        SnakemakeUtils.dump_tool_outputs(blastp, output)

rule bacmet_filter_blastp:
    """
    Filters the hits detected by blastp.
    """
    input:
        TSV = rules.bacmet_blastp.output.TSV,
        DB = rules.bacmet_pickle_db.output.DB
    output:
        TSV = Path(config['working_dir']) / 'bacmet' / 'hit_filtering' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'bacmet' / 'hit_filtering' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'bacmet' / 'hit_filtering',
        min_cov = 90,
        min_id = 75,
        cols = 'pident sseqid sseq slen qseqid qstart qend'
    run:
        from camel.app.io.tooliofile import ToolIOFile
        tsv_in = SnakemakeUtils.load_object(Path(input.TSV))[0].path
        data_in = pd.read_table(tsv_in, names=params.cols.split(' '))
        data_in['perc_covered'] = data_in.apply(lambda x: 100.0 * float(len(x['sseq'])) / x['slen'], axis=1)
        data_in_filt = data_in[data_in['pident'] > params.min_id].copy()
        data_in_filt = data_in_filt[data_in_filt['pident'] > params.min_cov].copy()
        data_in_filt['BacMet_ID'] = data_in_filt['sseqid'].apply(lambda x: x.split('|')[0])

        # Parse metadata
        dir_db = SnakemakeUtils.load_object(Path(input.DB))[0].path
        tsv_meta = next(dir_db.glob('*.tsv'))
        data_meta = pd.read_table(tsv_meta)
        data_out = pd.merge(data_in_filt, data_meta, on='BacMet_ID')
        path_out = Path(params.dir_, 'bacmet.tsv')
        data_out.to_csv(path_out, sep='\t', index=False)
        SnakemakeUtils.dump_object([ToolIOFile(path_out)], Path(output.TSV))
        SnakemakeUtils.dump_object(
            {'params': {'min_cov': params.min_cov, 'min_id': params.min_id}}, Path(output.INFORMS))

rule bacmet_report:
    """
    Creates the output report for BacMet. 
    """
    input:
        TSV = rules.bacmet_filter_blastp.output.TSV,
        INFORMS_blastp = rules.bacmet_blastp.output.INFORMS,
        INFORMS_filtering = rules.bacmet_filter_blastp.output.INFORMS,
        DB = rules.bacmet_pickle_db.output.DB
    output:
        HTML = Path(config['working_dir']) / 'bacmet' / 'report' / 'html.io'
    params:
        dir_ = Path(config['working_dir']) / 'bacmet' / 'report'
    run:
        from camel.app.tools.pipelines.klebsiella.bacmetreporter import BacMetReporter
        reporter = BacMetReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule bacmet_report_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'bacmet' / 'report' / 'html-empty.io'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('BacMet', Path(output.VAL_HTML))

rule bacmet_create_summary:
    """
    Creates a tabular summary output file for the BacMet assay.
    """
    input:
        TSV = rules.bacmet_filter_blastp.output.TSV
    output:
        TSV = Path(config['working_dir']) / 'bacmet' / 'summary_bacmet.tsv'
    run:
        path_tsv = SnakemakeUtils.load_object(Path(input.TSV))[0].path
        data_in = pd.read_table(path_tsv)
        with open(output.TSV, 'w') as handle:
            handle.write('\t'.join(['bacmet_genes', ', '.join(sorted(list(data_in['Gene_name'])))]))
            handle.write('\n')
