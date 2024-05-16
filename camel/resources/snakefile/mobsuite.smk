from pathlib import Path

from pandas.errors import EmptyDataError

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import mobsuite


rule mobsuite_mob_recon:
    """
    Runs the MOB-recon tool.
    """
    input:
        FASTA = Path(config['working_dir']) / mobsuite.INPUT_MOBSUITE_FASTA,
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

rule mobsuite_create_summary:
    """
    Creates a tabular summary output for MOB-suite.
    """
    input:
        TSV = rules.mobsuite_mob_recon.output.TSV,
        INFORMS = rules.mobsuite_mob_recon.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / 'mob_suite' / 'summary_mob_suite.tsv'
    run:
        import pandas as pd

        # Parse TSV output
        path_tsv = SnakemakeUtils.load_object(Path(input.TSV))[0].path
        try:
            data_mobsuite = pd.read_table(path_tsv)
            if 'primary_cluster_id' not in data_mobsuite.columns:
                raise EmptyDataError
            primary_cluster_ids = list(data_mobsuite['primary_cluster_id'])
        except EmptyDataError:
            logger.info(f'No plasmids detected by MOB-suite')
            primary_cluster_ids = []

        # Parse informs
        informs_in = SnakemakeUtils.load_object(Path(input.INFORMS))

        # Create summary output
        with open(output.TSV, 'w') as handle:
            # Cluster IDs
            handle.write('\t'.join([
                'mob_suite_primary_cluster_ids',
                ', '.join(primary_cluster_ids) if len(primary_cluster_ids) > 0 else '-']))
            handle.write('\n')

            # Contigs classified as plasmids
            handle.write('\t'.join([
                'mob_suite_predicted_plasmid_contigs',
                ', '.join(ctg for ctg, status in informs_in['contig_report'].items() if status is not None)]))
            handle.write('\n')

            # Tool info
            handle.write('\t'.join(['mob_suite_tool_version', informs_in['_name']]))
            handle.write('\n')

rule mobsuite_report_genomic_context:
    """
    Reports the genomic context for detected genes.
    """
    input:
        TSV = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'input' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'input' / 'informs.io',
        INFORMS_mob_recon = rules.mobsuite_mob_recon.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'html.io'
    params:
        dir_ = Path(config['working_dir']) / 'mob_suite' / 'genomic_context',
        detection_method = config['detection_method'],
        input_type = config.get('input_type', 'illumina')
    run:
        from camel.app.tools.mobsuite.genomiccontext import GenomicContext

        genomic_context = GenomicContext(Camel.get_instance())
        genomic_context.add_input_files({'TSV': SnakemakeUtils.load_object(Path(input.TSV))})
        genomic_context.add_input_informs({
            'mob_recon': SnakemakeUtils.load_object(Path(input.INFORMS_mob_recon)),
            'dbs': SnakemakeUtils.load_object(Path(input.INFORMS))
        })
        genomic_context.update_parameters(
            detection_method=str(params.detection_method), input_type=str(params.input_type))
        step = Step(str(rule), genomic_context, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(genomic_context, output)

rule mobsuite_report_genomic_context_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'html-empty.io'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Genomic context', Path(output.VAL_HTML), 2)
