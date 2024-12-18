from pathlib import Path

from camel.resources.snakefile import trimming_illumina, gene_detection, trimming, variant_calling, variant_filtering, \
    contamination_check_kraken, quality_checks, sequence_typing, downsampling, amrfinder, mobsuite, quast, confindr, \
    core, trimming_ont, assembly, human_read_scrubbing, resfinder4, bacmet
from camel.scripts.klebsiellapipeline.snakefile import kleborate

#######################
# Included Snakefiles #
#######################
include: core.SNAKEFILE_CORE
include: human_read_scrubbing.SNAKEFILE_SCRUBBING
include: downsampling.SNAKEFILE_DOWNSAMPLING
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: trimming_ont.SNAKEFILE_TRIMMING_ONT
include: assembly.SNAKEFILE_ASSEMBLY
include: quast.SNAKEFILE_QUAST
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: confindr.SNAKEFILE_CONFINDR
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: variant_calling.SNAKEFILE_VARIANT_CALLING
include: variant_filtering.SNAKEFILE_VARIANT_FILTERING
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING
include: amrfinder.SNAKEFILE_AMRFINDER
include: resfinder4.SNAKEFILE_RESFINDER4
include: kleborate.SNAKEFILE_KLEBORATE
include: mobsuite.SNAKEFILE_MOB_SUITE
include: bacmet.SNAKEFILE_BACMET


#########
# Rules #
#########
rule all:
    """
    Rule to generate the required output files.
    """
    input:
        HTML = config['output_report'],
        TSV = config['output_tabular']

rule report_command_section:
    """
    Creates a report section with the commands used in the pipeline. 
    """
    input:
        INFORMS_scrubbing = human_read_scrubbing.get_command_informs(config),
        INFORMS_downsampling = downsampling.get_command_informs(config),
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_assembly = assembly.get_command_informs(config),
        INFORMS_quast = Path(config['working_dir']) /quast.OUTPUT_QUAST_INFORMS,
        INFORMS_busco = Path(config['working_dir']) / quast.OUTPUT_BUSCO_INFORMS,
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        INFORMS_confindr = confindr.get_command_informs(config),
        # INFORMS_mapping = quality_checks.get_mapping_rate_informs(config),
        # INFORMS_depth = quality_checks.get_depth_informs(config),
        INFORMS_variant_calling_all = variant_calling.get_command_informs(config) if 'variant_calling' in config['analyses'] else [],
        INFORMS_variant_filtering_all = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_INFORMS_ALL if 'variant_calling' in config['analyses'] else [],
        INFORMS_amrfnder = Path(config['working_dir']) / amrfinder.OUTPUT_AMRFINDER_INFORMS if 'amrfinder' in config['analyses'] else [],
        INFORMS_resfinder4 = Path(config['working_dir']) / resfinder4.OUTPUT_RESFINDER4_INFORMS if 'resfinder4' in config['analyses'] else [],
        INFORMS_mob_suite = Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_INFORMS  if 'mob_suite' in config['analyses'] else [],
        INFORMS_kleborate = Path(config['working_dir']) / kleborate.OUTPUT_KLEBORATE_INFORMS if 'kleborate' in config['analyses'] else [],
        INFORMS_vfdb_core = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        IFNORMS_bacmet = Path(config['working_dir']) / bacmet.OUTPUT_BACMET_INFORMS if 'bacmet' in config['analyses'] else [],
        IFNORMS_prodigal = Path(config['working_dir']) / bacmet.OUTPUT_PRODIGAL_INFORMS if 'bacmet' in config['analyses'] else []
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.export_command_section(input, Path(output.HTML), Path(params.dir_))

rule report_combine_all:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        reports_scrubbing = human_read_scrubbing.get_reports(config),
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        report_quast = Path(config['working_dir']) / quast.OUTPUT_QUAST_REPORT,
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_confindr = confindr.get_report(config),
        report_adv_qc = Path(config['working_dir']) / str(quality_checks.OUTPUT_QUALITY_CHECKS_REPORT).format(
            input_type=config['input_type']),
        report_variant = variant_calling.get_reports(config) if 'variant_calling' in config['analyses'] else [],
        # Species identification
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        # AMR detection
        report_amrfinder = Path(config['working_dir']) / (amrfinder.OUTPUT_AMRFINDER_REPORT if 'amrfinder' in config['analyses'] else amrfinder.OUTPUT_AMRFINDER_REPORT_EMPTY),
        report_resfinder4 = Path(config['working_dir']) / (resfinder4.OUTPUT_RESFINDER4_REPORT if 'resfinder4' in config['analyses'] else resfinder4.OUTPUT_RESFINDER4_REPORT_EMPTY),
        # Virulence gene detection
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        # Plasmid characterization
        report_plasmidfinder = gene_detection.get_gene_detection_report('plasmidfinder', config),
        report_mob_suite = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_REPORT_EMPTY),
        report_genomic_context = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'html.io',
        # Kleborate
        report_kleborate = Path(config['working_dir']) / (kleborate.OUTPUT_KLEBORATE_REPORT if 'kleborate' in config['analyses'] else kleborate.OUTPUT_KLEBORATE_REPORT_EMPTY),
        # BacMet
        report_prodigal = Path(config['working_dir']) / (bacmet.OUTPUT_PRODIGAL_REPORT if 'bacmet' in config['analyses'] else bacmet.OUTPUT_PRODIGAL_REPORT_EMPTY),
        report_bacmet = Path(config['working_dir']) / (bacmet.OUTPUT_BACMET_REPORT if 'bacmet' in config['analyses'] else bacmet.OUTPUT_BACMET_REPORT_EMPTY),
        # Typing
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config),
        report_scgmlst = sequence_typing.get_sequence_typing_report('scgmlst', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        # Report
        report_citations = Path(config['working_dir'], core.OUTPUT_HTML_CITATIONS),
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        output_dir= config['output_dir'],
        pipeline_info = config['pipeline'],
        input_dict = config['input'],
        input_type = config['input_type'],
        detection_method = config['detection_method']
    run:
        import datetime
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        # Add the header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            sample_name=params.sample_name,
            date=datetime.datetime.now(),
            pipeline_version=params.pipeline_info['version'],
            input_files=ReportPipeline.format_input_string(params.input_dict),
            input_type=params.input_type,
            detection_method=params.detection_method
        ))

        # Set up the report content structure
        report_structure = []
        ReportPipeline.add_content_scrubbing(
            report_structure, params.input_type, input.reports_scrubbing)
        ReportPipeline.add_content_trim_basic_qc(
            report_structure,params.input_type,input.reports_downsampling, input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        ReportPipeline.add_content_contamination_check(
            report_structure,params.input_type,input.reports_contamination, input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))
        if 'variant_calling' in config['analyses']:
            report_structure.append(('Variant calling', 'variant', [Path(input.report_variant)]))
        report_structure.extend([
            ('Species identification', 'species', [Path(input.report_rmlst)]),
            ('AMR detection', 'amr', [Path(input.report_amrfinder), Path(input.report_resfinder4)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (input.report_vfdb_core,)]),
            ('Kleborate', 'kleborate', [Path(input.report_kleborate)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_plasmidfinder, input.report_mob_suite, input.report_genomic_context)]),
            ('Biocide and metal resistance', 'bacmet', [Path(input.report_prodigal), Path(input.report_bacmet)]),
            ('Sequence typing', 'st', [Path(x) for x in (
                input.report_mlst, input.report_scgmlst, input.report_cgmlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ])
        SnakePipelineUtils.add_report_content(report, report_structure)

rule summary_combine_all:
    """
    Combines the summary information of several steps into a single TSV file.
    """
    input:
        Path(config['working_dir'], core.OUTPUT_TSV_SUMMARY_INIT),
        human_read_scrubbing.get_summaries(config),
        downsampling.get_summaries(config),
        trimming.get_summaries(config),
        Path(config['working_dir']) / quast.OUTPUT_QUAST_SUMMARY,
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        contamination_check_kraken.get_summaries(config),
        confindr.get_summary(config),
        variant_calling.get_summaries(config) if 'variant_calling' in config['analyses'] else [],
        # AMR detection
        Path(config['working_dir']) / amrfinder.OUTPUT_AMRFINDER_SUMMARY if 'amrfinder' in config['analyses'] else [],
        Path(config['working_dir']) / resfinder4.OUTPUT_RESFINDER4_SUMMARY if 'resfinder4' in config['analyses'] else [],
        # Virulence detection
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        # Kleborate
        Path(config['working_dir']) / kleborate.OUTPUT_KLEBORATE_SUMMARY if 'kleborate' in config['analyses'] else [],
        # Plasmid characterization
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else [],
        Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_SUMMARY if 'mob_suite' in config['analyses'] else [],
        # BacMet
        Path(config['working_dir']) / bacmet.OUTPUT_BACMET_SUMMARY if 'bacmet' in config['analyses'] else [],
        # Sequence typing
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='rmlst') if 'rmlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='scgmlst') if 'scgmlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else []
    output:
        TSV = config['output_tabular']
    run:
        with open(output.TSV, 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())

rule link_genomic_context:
    """
    Links the input databases to the genomic context assay.
    """
    input:
        # AMR
        TSV_amrfinder = Path(config['working_dir']) / 'amrfinder' / 'tsv.io' if 'amrfinder' in config['analyses'] else [],
        # Virulence
        TSV_gd_vfdb = Path(config['working_dir']) / 'gene_detection' / 'vfdb_core' / 'metadata' / 'tsv.io' if 'vfdb_core' in config['analyses'] else [],
        INFORMS_gd_vfdb = Path(config['working_dir']) / 'gene_detection' / 'vfdb_core' / 'db_manager' / 'informs.io' if 'vfdb_core' in config['analyses'] else [],
        # BacMet
        TSV_bacmet = Path(config['working_dir']) / 'bacmet' / 'hit_filtering' / 'tsv.io' if 'bacmet' in config['analyses'] else []
    output:
        TSV = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'input' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'input' / 'informs.io'
    run:
        mobsuite.collect_genomic_context_input(input, Path(output.TSV), Path(output.INFORMS))
