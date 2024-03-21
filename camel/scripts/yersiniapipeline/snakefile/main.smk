from pathlib import Path

from camel.resources.snakefile import trimming, trimming_illumina, quality_checks, \
    contamination_check_kraken, gene_detection, sequence_typing, downsampling, confindr, quast, core, trimming_ont, mobsuite, resfinder4, amrfinder, \
    assembly
from camel.scripts.yersiniapipeline.snakefile import species_determination

#######################
# Included Snakefiles #
#######################

include: core.SNAKEFILE_CORE
include: downsampling.SNAKEFILE_DOWNSAMPLING
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: trimming_ont.SNAKEFILE_TRIMMING_ONT
include: assembly.SNAKEFILE_ASSEMBLY
include: quast.SNAKEFILE_QUAST
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: confindr.SNAKEFILE_CONFINDR
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: amrfinder.SNAKEFILE_AMRFINDER
include: resfinder4.SNAKEFILE_RESFINDER4
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: mobsuite.SNAKEFILE_MOB_SUITE
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING
include: species_determination.SNAKEFILE_SPECIES_DETERMINATION

#########
# Rules #
#########

rule all:
    """
    This rule ensures that the required output files are generated.
    """
    input:
        config['output_report'],
        config['output_tabular']

rule report_create_command_section:
    """
    Creates the report section containing the tool commands.
    """
    input:
        INFORMS_downsampling = downsampling.get_command_informs(config),
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_assembly = assembly.get_command_informs(config),
        INFORMS_quast = Path(config['working_dir']) / quast.OUTPUT_QUAST_INFORMS,
        INFORMS_busco = Path(config['working_dir']) / quast.OUTPUT_BUSCO_INFORMS,
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        INFORMS_confindr = confindr.get_command_informs(config),
        INFORMS_assembly_map = assembly.get_qc_informs(config,config['input_type']),
        INFORMS_amrfinder = Path(config['working_dir']) / amrfinder.OUTPUT_AMRFINDER_INFORMS if 'amrfinder' in config['analyses'] else [],
        INFORMS_resfinder4=Path(config['working_dir']) / resfinder4.OUTPUT_RESFINDER4_INFORMS if 'resfinder4' in config[
            'analyses'] else [],
        INFORMS_mob_suite = Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_INFORMS if 'mob_suite' in config['analyses'] else [],
        INFORMS_vfdb_core = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        INFORMS_mlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        INFORMS_mlst_mcnally = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='mlst_mcnally') if 'mlst_mcnally' in config['analyses'] else [],
        INFORMS_amr = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='resistance_genes') if 'resistance_genes' in config['analyses'] else [],
        INFORMS_cgmlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        INFORMS_cgmlst_ye = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='cgmlst_ye') if 'cgmlst_ye' in config['analyses'] else [],
        INFORMS_cgmlst_yp = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='cgmlst_yp') if 'cgmlst_yp' in config['analyses'] else [],
        INFORMS_cgmlst_enterobase = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='cgmlst_enterobase') if 'cgmlst_enterobase' in config['analyses'] else [],
        INFORMS_rmlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='rmlst') if 'rmlst' in config['analyses'] else []
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.export_command_section(input, Path(output.HTML), Path(params.dir_))

rule combine_reports:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        report_quast = Path(config['working_dir']) / quast.OUTPUT_QUAST_REPORT,
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_confindr = confindr.get_report(config),
        report_adv_qc = Path(config['working_dir']) / str(quality_checks.OUTPUT_QUALITY_CHECKS_REPORT).format(
            input_type=config['input_type']),
        # AMR detection
        report_amrfinder = Path(config['working_dir']) / (amrfinder.OUTPUT_AMRFINDER_REPORT if 'amrfinder' in config['analyses'] else amrfinder.OUTPUT_AMRFINDER_REPORT_EMPTY),
        report_resfinder4 = Path(config['working_dir']) / (resfinder4.OUTPUT_RESFINDER4_REPORT if 'resfinder4' in config['analyses'] else resfinder4.OUTPUT_RESFINDER4_REPORT_EMPTY),
        # Virulence gene detection
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core',config),
        # Genomic context investigation
        report_mob_suite = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_REPORT_EMPTY),
        report_genomic_context = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_CONTEXT_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_CONTEXT_REPORT_EMPTY),
        # Sequence typing
        report_mlst = sequence_typing.get_sequence_typing_report('mlst',config),
        report_mlst_mcnally = sequence_typing.get_sequence_typing_report('mlst_mcnally',config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst',config),
        report_cgmlst_ye = sequence_typing.get_sequence_typing_report('cgmlst_ye',config),
        report_cgmlst_yp = sequence_typing.get_sequence_typing_report('cgmlst_yp',config),
        report_cgmlst_enterobase = sequence_typing.get_sequence_typing_report('cgmlst_enterobase',config),
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        # Species determination
        report_species = Path(config['working_dir']) / (species_determination.OUTPUT_SPECIES_DETERMINATION_REPORT if 'species' in config['analyses'] else species_determination.OUTPUT_SPECIES_DETERMINATION_REPORT_EMPTY),
        # Report
        report_citations = Path(config['working_dir'], core.OUTPUT_HTML_CITATIONS),
        report_commands = rules.report_create_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        input_dict = config['input'],
        input_type = config['input_type'],
        citation_keys = config['citations'],
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
            key_citation=params.citation_keys['main'],
            detection_method=params.detection_method
        ))

        # Set up the report content structure
        report_structure = []
        ReportPipeline.add_content_trim_basic_qc(
            report_structure, params.input_type, input.reports_downsampling, input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        ReportPipeline.add_content_contamination_check(
            report_structure, params.input_type, input.reports_contamination, input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))
        report_structure.extend([
            ('AMR detection', 'amr', [Path(x) for x in (
                input.report_amrfinder, input.report_resfinder4)]),
            ('Virulence detection', 'virulence', [Path(input.report_vfdb_core)]),
            ('Genomic context', 'mob_suite', [Path(x) for x in(
                input.report_mob_suite, input.report_genomic_context)]),
            ('Sequence typing', 'st', [Path(x) for x in (
                input.report_mlst, input.report_mlst_mcnally, input.report_cgmlst, input.report_cgmlst_ye, input.report_cgmlst_yp, input.report_cgmlst_enterobase, input.report_rmlst)]),
            ('Species determination', 'species', [Path(input.report_species)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ])
        SnakePipelineUtils.add_report_content(report, report_structure)


rule combine_summary_files:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        Path(config['working_dir'], core.OUTPUT_TSV_SUMMARY_INIT),
        downsampling.get_summaries(config),
        trimming.get_summaries(config),
        Path(config['working_dir']) / quast.OUTPUT_QUAST_SUMMARY,
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        contamination_check_kraken.get_summaries(config),
        confindr.get_summary(config),
        Path(config['working_dir']) / amrfinder.OUTPUT_AMRFINDER_SUMMARY if 'amrfinder' in config['analyses'] else [],
        Path(config['working_dir']) / resfinder4.OUTPUT_RESFINDER4_SUMMARY if 'resfinder4' in config['analyses'] else [],
        Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_SUMMARY if 'mob_suite' in config['analyses'] else [],
        # Virulence detection
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        # Sequence typing
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst_mcnally') if 'mlst_mcnally' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst_ye') if 'cgmlst_ye' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst_yp') if 'cgmlst_yp' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst_enterobase') if 'cgmlst_enterobase' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='rmlst') if 'rmlst' in config['analyses'] else [],
        Path(config['working_dir']) / species_determination.OUTPUT_SPECIES_DETERMINATION_SUMMARY if 'species' in config['analyses'] else []
    output:
        TSV = config.get('output_tabular')
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
        INFORMS_gd_vfdb = Path(config['working_dir']) / 'gene_detection' / 'vfdb_core' / 'db_manager' / 'informs.io' if 'vfdb_core' in config['analyses'] else []
    output:
        TSV = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'input' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'input' / 'informs.io'
    run:
        mobsuite.collect_genomic_context_input(input, Path(output.TSV), Path(output.INFORMS))

rule link_species:
    """
    Links the cgST profile matches to the species determination tool.
    """
    input:
        TSV = Path(config['working_dir']) / 'typing' / 'cgmlst' / 'detect_st' / 'tsv.io'
    output:
        TSV_profile_matches = Path(config['working_dir']) / 'species_determination' / 'input' / 'tsv_profile_matches.io',
        TSV_taxonomic = Path(config['working_dir']) / 'species_determination' / 'input' / 'tsv_taxonomic.io'
    params:
        TSV = Path(config['species_determination']['tsv'])
    run:
        import shutil
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        from camel.app.io.tooliofile import ToolIOFile
        shutil.copyfile(input.TSV, output.TSV_profile_matches)
        SnakemakeUtils.dump_object([ToolIOFile(params.TSV)], Path(output.TSV_taxonomic))
