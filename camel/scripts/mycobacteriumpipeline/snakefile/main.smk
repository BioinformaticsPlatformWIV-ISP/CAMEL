from pathlib import Path

from camel.resources.snakefile import trimming_illumina, contamination_check_kraken, quality_checks, variant_calling, \
    variant_filtering, gene_detection, sequence_typing, trimming, downsampling, confindr, quast, core, assembly, \
    human_read_scrubbing, read_simulation
from camel.scripts.mycobacteriumpipeline.snakefile import csb_rd, snpit, hsp65, spoligotyping, snplineage, assay51snp, \
    amrdetection

#######################
# Included Snakefiles #
#######################
include: core.SNAKEFILE_CORE
include: human_read_scrubbing.SNAKEFILE_SCRUBBING
include: read_simulation.SNAKEFILE_READ_SIMULATION
include: downsampling.SNAKEFILE_DOWNSAMPLING
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: assembly.SNAKEFILE_ASSEMBLY
include: confindr.SNAKEFILE_CONFINDR
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: quast.SNAKEFILE_QUAST
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING
include: variant_calling.SNAKEFILE_VARIANT_CALLING
include: variant_filtering.SNAKEFILE_VARIANT_FILTERING
include: csb_rd.SNAKEFILE_CSB_RD
include: snpit.SNAKEFILE_SNPIT
include: hsp65.SNAKEFILE_HSP65
include: assay51snp.SNAKEFILE_51SNP
include: spoligotyping.SNAKEFILE_SPOLIGOTYPING
include: snplineage.SNAKEFILE_SNP_LINEAGE
include: amrdetection.SNAKEFILE_AMR

#########
# Rules #
#########
rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        config['output_report'],
        config['output_tabular']

rule report_command_section:
    """
    Creates a HTML report section with the main commands.
    """
    input:
        INFORMS_scrubbing = human_read_scrubbing.get_command_informs(config),
        INFORMS_downsampling = downsampling.get_command_informs(config),
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_assembly = assembly.get_command_informs(config),
        INFORMS_quast = Path(config['working_dir']) / quast.OUTPUT_QUAST_INFORMS,
        INFORMS_busco = Path(config['working_dir']) / quast.OUTPUT_BUSCO_INFORMS,
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        INFORMS_confindr = confindr.get_command_informs(config),
        INFORMS_assembly_map = assembly.get_qc_informs(config, config['input_type'], mode='ref'),
        INFORMS_variant_calling_all = variant_calling.get_command_informs(config),
        INFORMS_variant_filtering_all = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_INFORMS_ALL,
        INFORMS_snpit = Path(config['working_dir']) / snpit.OUTPUT_SNPIT_INFORMS if 'snpit' in config['analyses'] else [],
        INFORMS_16s = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='ncbi_16s') if 'ncbi_16s' in config['analyses'] else [],
        INFORMS_csb_rd = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='csb_rd') if 'csb_rd' in config['analyses'] else [],
        INFORMS_hsp65 = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='hsp65') if 'hsp65' in config['analyses'] else [],
        INFORMS_spoligo = Path(config['working_dir']) / spoligotyping.OUTPUT_SPOLIGOTYPING_INFORMS if 'spoligotyping' in config['analyses'] else [],
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.export_command_section(input, Path(output.HTML), Path(params.dir_))

rule report_combine_all:
    """
    Combines the HTML report sections into a single report.
    """
    input:
        reports_scrubbing = human_read_scrubbing.get_reports(config),
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        report_quast = Path(config['working_dir']) / quast.OUTPUT_QUAST_REPORT,
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_confindr = confindr.get_report(config),
        report_adv_qc = Path(config['working_dir']) / str(quality_checks.OUTPUT_QUALITY_CHECKS_REPORT).format(input_type=config['input_type']),
        report_variant = variant_calling.get_reports(config),
        # Species identification
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        report_ncbi_16s = gene_detection.get_gene_detection_report('ncbi_16s', config),
        report_51snp = Path(config['working_dir']) / (assay51snp.OUTPUT_51SNP_REPORT if '51snp' in config['analyses'] else assay51snp.OUTPUT_51SNP_REPORT_EMPTY),
        report_csb_rd = Path(config['working_dir']) / (csb_rd.OUTPUT_CSB_RD_REPORT if 'csb_rd' in config['analyses'] else csb_rd.OUTPUT_CSB_RD_REPORT_EMPTY),
        report_hsp65 = Path(config['working_dir']) / (hsp65.OUTPUT_HSP65_REPORT if 'hsp65' in config['analyses'] else hsp65.OUTPUT_HSP65_REPORT_EMPTY),
        report_snpit = Path(config['working_dir']) / (snpit.OUTPUT_SNPIT_REPORT if 'snpit' in config['analyses'] else snpit.OUTPUT_SNPIT_REPORT_EMPTY),
        # Spoligotyping & lineage determination
        report_spoligo = Path(config['working_dir']) / (spoligotyping.OUTPUT_SPOLIGOTYPING_REPORT if 'spoligotyping' in config['analyses'] else spoligotyping.OUTPUT_SPOLIGOTYPING_REPORT_EMPTY),
        report_snp_lineage = Path(config['working_dir']) / (snplineage.OUTPUT_SNP_LINEAGE_REPORT if 'snp_lineage' in config['analyses'] else snplineage.OUTPUT_SNP_LINEAGE_REPORT_EMPTY),
        # AMR
        report_amr = Path(config['working_dir']) / (amrdetection.OUTPUT_AMR_REPORT if 'amr' in config['analyses'] else amrdetection.OUTPUT_AMR_REPORT_EMPTY),
        report_amr_genes = Path(config['working_dir']) / 'amr' / 'cds' / 'html.io',
        # Typing
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        # Report
        report_citations = Path(config['working_dir'], core.OUTPUT_HTML_CITATIONS),
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        input_dict = config['input'],
        input_type = config['input_type'],
        citation_keys = config['citations']
    run:
        import datetime
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        # Add the header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        section = SnakePipelineUtils.create_input_section(
            sample_name=params.sample_name,
            date=datetime.datetime.now(),
            pipeline_version=params.pipeline_info['version'],
            input_files=ReportPipeline.format_input_string(params.input_dict),
            input_type=params.input_type,
            key_citation=params.citation_keys['main']
        )
        if params.input_type == 'fasta':
            section.add_warning_message(
                'SNP-based assays are run on simulated reads from the assembled contigs, which may differ from the '
                'original reads.')
        report.add_html_object(section)

        # Add report content
        report_structure = []
        ReportPipeline.add_content_scrubbing(
            report_structure, params.input_type, input.reports_scrubbing)
        ReportPipeline.add_content_trim_basic_qc(
            report_structure,params.input_type,input.reports_downsampling,input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        ReportPipeline.add_content_contamination_check(
            report_structure,params.input_type,input.reports_contamination,input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))

        # Add output sections
        report_structure.extend([
            ('Variant calling', 'variant', [Path(input.report_variant)]),
            ('Species identification', 'identification', [
                Path(input.report_rmlst), Path(input.report_ncbi_16s), Path(input.report_snpit),
                Path(input.report_csb_rd), Path(input.report_hsp65), Path(input.report_51snp)]),
            ('Spoligotyping and lineage', 'spoligotyping', [
                Path(input.report_spoligo), Path(input.report_snp_lineage)]),
            ('AMR detection', 'amr', [Path(input.report_amr), Path(input.report_amr_genes)]),
            ('Sequence typing', 'typing', [Path(input.report_mlst), Path(input.report_cgmlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ])
        SnakePipelineUtils.add_report_content(report, report_structure)

rule summary_combine_all:
    """
    Combines the summary output files of the different assays into a single summary file.
    """
    input:
        Path(config['working_dir'], core.OUTPUT_TSV_SUMMARY_INIT),
        human_read_scrubbing.get_summaries(config),
        downsampling.get_summaries(config),
        trimming.get_summaries(config),
        Path(config['working_dir']) / quast.OUTPUT_QUAST_SUMMARY,
        contamination_check_kraken.get_summaries(config),
        confindr.get_summary(config),
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        variant_calling.get_summaries(config),
        Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_SUMMARY,
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='ncbi_16s') if 'ncbi_16s' in config['analyses'] else [],
        Path(config['working_dir']) / csb_rd.OUTPUT_CSB_RD_SUMMARY if 'csb_rd' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='hsp65') if 'hsp65' in config['analyses'] else [],
        Path(config['working_dir']) / assay51snp.OUTPUT_51SNP_SUMMARY if '51snp' in config['analyses'] else [],
        Path(config['working_dir']) / snpit.OUTPUT_SNPIT_SUMMARY if 'snpit' in config['analyses'] else [],
        Path(config['working_dir']) / spoligotyping.OUTPUT_SPOLIGOTYPING_SUMMARY if 'spoligotyping' in config['analyses'] else [],
        Path(config['working_dir']) / snplineage.OUTPUT_SNP_LINEAGE_SUMMARY if 'snp_lineage' in config['analyses'] else [],
        Path(config['working_dir']) / amrdetection.OUTPUT_AMR_SUMMARY if 'amr' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='rmlst') if 'rmlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else []
    output:
        TSV = config['output_tabular']
    run:
        with open(output.TSV, 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
