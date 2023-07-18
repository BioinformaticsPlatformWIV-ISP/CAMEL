from pathlib import Path

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
        path_tsv = SnakemakeUtils.load_object(Path(input.TSV))[0].path
        data_mobsuite = pd.read_table(path_tsv)

        informs_in = SnakemakeUtils.load_object(Path(input.INFORMS))
        with open(output.TSV, 'w') as handle:
            # Primary cluster id
            try:
                handle.write('\t'.join([
                    'mob_suite_primary_cluster_ids', ', '.join(list(data_mobsuite['primary_cluster_id']))]))
                handle.write('\n')

                # Contigs classified as plasmids
                handle.write('\t'.join([
                    'mob_suite_predicted_plasmid_contigs',
                    ', '.join(ctg for ctg, status in informs_in['contig_report'].items() if status is not None)
                ]))
                handle.write('\n')
            except KeyError:
                handle.write('No plasmids found.')
                handle.write('\n')

rule mobsuite_report_genomic_context:
    """
    Reports the genomic context for detected genes.
    """
    input:
        TSV_amrfinder = Path(config['working_dir']) / 'amrfinder' / 'tsv.io' if 'amrfinder' in config['analyses'] else [],
        TSV_bacmet = Path(config['working_dir']) / 'bacmet' / 'hit_filtering' / 'tsv.io' if 'bacmet' in config['analyses'] else [],
        TSV_vfdb = Path(config['working_dir']) / 'gene_detection' / 'vfdb_core' / 'metadata' / 'tsv.io' if 'vfdb_core' in config['analyses'] else [],
        INFORMS_mob_recon = rules.mobsuite_mob_recon.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'html.io'
    params:
        dir_ = Path(config['working_dir']) / 'mob_suite' / 'genomic_context',
        detection_method = config['detection_method'],
        read_type = config.get('read_type', 'illumina')
    run:
        from camel.app.tools.pipelines.klebsiella.genomiccontext import GenomicContext
        genomic_context = GenomicContext(Camel.get_instance())
        informs = {
            'TSV_amrfinder': {'key': 'amrfinder', 'title': 'AMRFinder', 'contig': 'Contig id', 'gene': 'Gene symbol'},
            'TSV_bacmet': {'key': 'bacmet', 'title': 'BacMet', 'contig': 'qseqid', 'gene': 'Gene_name'},
            'TSV_vfdb': {'key': 'vfdb', 'title': 'VFDB core', 'contig': 'Contig', 'gene': 'Gene'}
        }
        db_informs_to_add = []
        for k, v in input.items():
            if not v:
                continue
            if 'TSV' in k:
                db_informs_to_add.append(informs[k])
                SnakemakeUtils.add_pickle_input(genomic_context, k, Path(v))
        genomic_context.add_input_informs({
            'mob_recon': SnakemakeUtils.load_object(Path(input.INFORMS_mob_recon)),
            'dbs': db_informs_to_add})
        genomic_context.update_parameters(
            detection_method=str(params.detection_method), read_type=str(params.read_type))
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
        SnakePipelineUtils.create_empty_report_section('Genomic context', Path(output.VAL_HTML), 3)
