from pathlib import Path
from camel.resources.snakefile import (trimming, trimming_illumina, trimming_ont, quality_checks,
contamination_check_kraken, gene_detection, sequence_typing, downsampling, quast, confindr, core, assembly, resfinder4,
mobsuite, abritamr, mykrobe, human_read_scrubbing)
from camel.scripts.salmonellapipeline.snakefile import spifinder, serotyping_salmonella

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
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: resfinder4.SNAKEFILE_RESFINDER4
include: mobsuite.SNAKEFILE_MOB_SUITE
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING
include: serotyping_salmonella.SNAKEFILE_SEROTYPE
include: mykrobe.SNAKEFILE_MYKROBE
include: spifinder.SNAKEFILE_SPIFINDER
include: abritamr.SNAKEFILE_ABRITAMR

#########
# Rules #
#########
rule all:
    input:
        HTML = config['output_report'],
        TSV = config['output_tabular'],
        #JSON = config['output_json']

rule report_command_section:
    input:
        INFORMS_human_read_scrubbing = human_read_scrubbing.get_command_informs(config),
        INFORMS_downsampling = downsampling.get_command_informs(config),
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_assembly = assembly.get_command_informs(config),
        INFORMS_quast = Path(config['working_dir']) / quast.OUTPUT_QUAST_INFORMS,
        INFORMS_busco = Path(config['working_dir']) / quast.OUTPUT_BUSCO_INFORMS,
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        INFORMS_confindr = confindr.get_command_informs(config),
        INFORMS_assembly_map = assembly.get_qc_informs(config,config['input_type']),
        #INFORMS_variant_calling_all = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_INFORMS_ALL,
        #INFORMS_variant_filtering_all = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_INFORMS_ALL,
        INFORMS_serotyping = serotyping_salmonella.get_command_informs(config),
        INFORMS_mykrobe = Path(config['working_dir']) / mykrobe.OUTPUT_MYKROBE_INFORMS if 'mykrobe' in config['analyses'] else [],
        INFORMS_abritamr_run =  abritamr.get_command_informs(config),
        INFORMS_resfinder4 = Path(config['working_dir']) / resfinder4.OUTPUT_RESFINDER4_INFORMS if 'resfinder4' in config['analyses'] else [],
        INFORMS_spifinder = spifinder.get_command_informs(config),
        INFORMS_vfdb_core = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        INFORMS_mob_suite = Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_INFORMS if 'mob_suite' in config['analyses'] else [],
        INFORMS_mlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        INFORMS_cgmlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        INFORMS_rmlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='rmlst') if 'rmlst' in config['analyses'] else []
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.export_command_section(input,Path(output.HTML),Path(params.dir_))


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
        report_adv_qc = Path(config['working_dir']) / str(quality_checks.OUTPUT_QUALITY_CHECKS_REPORT).format(input_type=config['input_type']),
        #report_variant = Path(config['working_dir']) /variant_calling.OUTPUT_VARIANT_CALLING_REPORT,
        # Gene detection
        report_resfinder4 = Path(config['working_dir']) / (resfinder4.OUTPUT_RESFINDER4_REPORT if 'resfinder4' in config['analyses'] else resfinder4.OUTPUT_RESFINDER4_REPORT_EMPTY),
        reports_spifinder = spifinder.get_reports(config),
        report_mykrobe = Path(config['working_dir']) / (mykrobe.OUTPUT_MYKROBE_REPORT if 'mykrobe' in config['analyses'] else mykrobe.OUTPUT_MYKROBE_REPORT_EMPTY),
        report_abritamr = Path(config['working_dir']) / (abritamr.OUTPUT_ABRITAMR_REPORT if 'abritamr' in config['analyses'] else abritamr.OUTPUT_ABRITAMR_REPORT_EMPTY),
        reports_serotyping = serotyping_salmonella.get_reports(config),
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        report_mob_suite = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_REPORT_EMPTY),
        report_genomic_context = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_CONTEXT_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_CONTEXT_REPORT_EMPTY),
        # Typing
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
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

        # Add report content
        report_structure = []
        ReportPipeline.add_content_scrubbing(
            report_structure, params.input_type, input.reports_scrubbing)
        ReportPipeline.add_content_trim_basic_qc(
            report_structure,params.input_type,input.reports_downsampling,input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        ReportPipeline.add_content_contamination_check(
            report_structure,params.input_type,input.reports_contamination,input.report_confindr)
        report_structure.extend([
            ('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]),
            ('Species identification', 'species', [Path(input.report_rmlst)]),
            ])
        ReportPipeline.add_content_serotyping_salmonella(
            report_structure,params.input_type,input.reports_serotyping)
        report_structure.extend([
            ('Lineage identification', 'mykrobe', [Path(input.report_mykrobe)]),
            ('AMR detection', 'amr', [Path(x) for x in (
                input.report_abritamr, input.report_resfinder4)])
             ])
        ReportPipeline.add_content_spifinder(
            report_structure,params.input_type,input.reports_spifinder)
        report_structure.extend([
            ('Virulence characterization', 'virulence', [Path(input.report_vfdb_core)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_mob_suite, input.report_genomic_context)]),
            ('Sequence typing', 'st', [Path(x) for x in (input.report_mlst, input.report_cgmlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ])
        SnakePipelineUtils.add_report_content(report, report_structure)

rule summary_combine_all:
    """
    Combines the summary information of several steps into a single TSV file.
    """
    input:
        Path(config['working_dir'],core.OUTPUT_TSV_SUMMARY_INIT),
        human_read_scrubbing.get_summaries(config),
        downsampling.get_summaries(config),
        trimming.get_summaries(config),
        Path(config['working_dir']) / quast.OUTPUT_QUAST_SUMMARY,
        contamination_check_kraken.get_summaries(config),
        confindr.get_summary(config),
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        #Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_SUMMARY,
        #Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_SUMMARY,
        serotyping_salmonella.get_summaries(config),
        Path(config['working_dir']) / mykrobe.OUTPUT_MYKROBE_SUMMARY if 'mykrobe' in config['analyses'] else [],
        Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_SUMMARY if 'abritamr' in config['analyses'] else [],
        Path(config['working_dir']) / resfinder4.OUTPUT_RESFINDER4_SUMMARY if 'resfinder4' in config['analyses'] else [],
        spifinder.get_summaries(config),
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_SUMMARY if 'mob_suite' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='rmlst') if 'rmlst' in config['analyses'] else []
    output:
        config.get('output_tabular')
    run:
        with open(output[0], 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())


rule link_genomic_context:
    """
    Links the input databases to the genomic context assay.
    """
    input:
        # AMR
        TSV_amrfinder = Path(config['working_dir']) / 'abritamr' / 'abritamr_output_amrfinder.io',
        # Virulence
        TSV_gd_vfdb = Path(config['working_dir']) / 'gene_detection' / 'vfdb_core' / 'metadata' / 'tsv.io' if 'vfdb_core' in config['analyses'] else [],
        INFORMS_gd_vfdb = Path(config['working_dir']) / 'gene_detection' / 'vfdb_core' / 'db_manager' / 'informs.io' if 'vfdb_core' in config['analyses'] else []
    output:
        TSV = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'input' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'input' / 'informs.io'
    run:
        mobsuite.collect_genomic_context_input(input, Path(output.TSV), Path(output.INFORMS))


# rule summary_json_combine_all:
#     """
#     Combines the summary information of several steps into a single JSON file.
#     """
#     input:
#         Path(config['working_dir'], core.OUTPUT_JSON_SUMMARY_INIT),
#         human_read_scrubbing.get_jsons(config),
#         downsampling.get_jsons(config),
#         trimming.get_jsons(config),
#         Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_SUMMARY_JSON,
#         Path(config['working_dir']) / quast.OUTPUT_QUAST_SUMMARY_JSON,
#         contamination_check_kraken.get_jsons(config),
#         confindr.get_json(config),
#         Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY_JSON,
#         #Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_SUMMARY_JSON if 'fasta' not in config['input'] else [],
#         #Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_SUMMARY_JSON if 'fasta' not in config['input'] else [],
#         Path(config['working_dir']) / resfinder4.OUTPUT_RESFINDER4_SUMMARY_JSON if 'resfinder4' in config['analyses'] else [],
#         Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_SUMMARY_JSON if 'mob_suite' in config['analyses'] else [],
#         Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY_JSON).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
#         Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY_JSON).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
#         spifinder.get_jsons(config),
#         Path(config['working_dir']) / mykrobe.OUTPUT_MYKROBE_SUMMARY_JSON if 'mykrobe' in config['analyses'] else [],
#         Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_SUMMARY_JSON if 'abritamr' in config['analyses'] else [],
#         serotyping_salmonella.get_jsons(config),
#         Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY_JSON).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
#         Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY_JSON).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
#         Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY_JSON).format(scheme='rmlst') if 'rmlst' in config['analyses'] else []
#     output:
#         JSON = config['output_json']
#     run:
#         import json
#         with open(output.JSON, 'w') as handle_out:
#             with open(input[0], 'r') as summary_metadata:
#                 mega_dict = json.load(summary_metadata)
#                 for summary_input in input[1:]:
#                     with open(summary_input, 'r') as assay_results_dictionary:
#                         mega_dict.update(json.load(assay_results_dictionary))
#                 handle_out.write(json.dumps(mega_dict))
