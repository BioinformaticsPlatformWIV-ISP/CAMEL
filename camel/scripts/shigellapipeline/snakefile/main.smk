from pathlib import Path

from camel.resources.snakefile import trimming, trimming_illumina, \
    quality_checks, contamination_check_kraken, gene_detection, sequence_typing, \
    downsampling, confindr, quast, core, assembly, amrfinder, resfinder4, mobsuite, mykrobe, human_read_scrubbing
from camel.scripts.shigellapipeline.snakefile import shigeifinder, shigatyper

#######################
# Included Snakefiles #
#######################
include: core.SNAKEFILE_CORE
include: human_read_scrubbing.SNAKEFILE_SCRUBBING
include: downsampling.SNAKEFILE_DOWNSAMPLING
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: assembly.SNAKEFILE_ASSEMBLY
include: quast.SNAKEFILE_QUAST
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: confindr.SNAKEFILE_CONFINDR
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: amrfinder.SNAKEFILE_AMRFINDER
include: resfinder4.SNAKEFILE_RESFINDER4
include: mobsuite.SNAKEFILE_MOB_SUITE
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING
include: shigeifinder.SNAKEFILE_SHIGEIFINDER
include: shigatyper.SNAKEFILE_SHIGATYPER
include: mykrobe.SNAKEFILE_MYKROBE


#########
# Rules #
#########
rule all:
    """
    This rules ensures that the required output files are generated.
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
        INFORMS_busco = Path(config['working_dir']) /quast.OUTPUT_BUSCO_INFORMS,
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        INFORMS_confindr = confindr.get_command_informs(config),
        INFORMS_assembly_map = assembly.get_qc_informs(config, config['input_type']),
        INFORMS_shigeifinder = Path(config['working_dir']) / shigeifinder.OUTPUT_SHIGEIFINDER_INFORMS if 'shigeifinder' in config['analyses'] else [],
        INFORMS_shigatyper = Path(config['working_dir']) / shigatyper.OUTPUT_SHIGATYPER_INFORMS if 'shigatyper' in config['analyses'] else[],
        INFORMS_mykrobe = Path(config['working_dir']) / mykrobe.OUTPUT_MYKROBE_INFORMS if 'mykrobe' in config['analyses'] else[],
        INFORMS_amrfinder = Path(config['working_dir']) / amrfinder.OUTPUT_AMRFINDER_INFORMS if 'amrfinder' in config['analyses'] else [],
        INFORMS_resfinder4 = Path(config['working_dir']) / resfinder4.OUTPUT_RESFINDER4_INFORMS if 'resfinder4' in config['analyses'] else [],
        INFORMS_virulence = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        INFORMS_virulence_shiga = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='virulencefinder_shiga') if 'virulencefinder' in config['analyses'] else [],
        INFORMS_serotype_o = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='serotype_o') if 'serotype_o' in config['analyses'] else[],
        INFORMS_mob_suite = Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_INFORMS if 'mob_suite' in config['analyses'] else[],
        INFORMS_mlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='mlst') if 'mlst' in config['analyses'] else[],
        INFORMS_rmlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='rmlst') if 'rmlst' in config['analyses'] else[]
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.export_command_section(input,Path(output.HTML),Path(params.dir_))

rule combine_reports:
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
        report_adv_qc = Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT,
        # Species identification
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst',config),
        # Gene detection
        report_amrfinder = Path(config['working_dir']) / (amrfinder.OUTPUT_AMRFINDER_REPORT if 'amrfinder' in config['analyses'] else amrfinder.OUTPUT_AMRFINDER_REPORT_EMPTY),
        report_resfinder4 = Path(config['working_dir']) / (resfinder4.OUTPUT_RESFINDER4_REPORT if 'resfinder4' in config['analyses'] else resfinder4.OUTPUT_RESFINDER4_REPORT_EMPTY),
        report_virulence = gene_detection.get_gene_detection_report('virulencefinder', config),
        report_virulence_shiga = gene_detection.get_gene_detection_report('virulencefinder_shiga', config, 'virulencefinder'),
        report_serotype_o_type= gene_detection.get_gene_detection_report('serotype_o',config,'serotype_o'),
        # Plasmid characterization
        report_mob_suite = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_REPORT_EMPTY),
        report_genomic_context = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_CONTEXT_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_CONTEXT_REPORT_EMPTY),
        # Shigella serotyping
        report_shigeifinder = Path(config['working_dir']) / (shigeifinder.OUTPUT_SHIGEIFINDER_REPORT if 'shigeifinder' in config['analyses'] else shigeifinder.OUTPUT_SHIGEIFINDER_REPORT_EMPTY),
        report_shigatyper = Path(config['working_dir']) / (shigatyper.OUTPUT_SHIGATYPER_REPORT if 'shigatyper' in config['analyses'] else shigatyper.OUTPUT_SHIGATYPER_REPORT_EMPTY),
        report_mykrobe = Path(config['working_dir']) / (mykrobe.OUTPUT_MYKROBE_REPORT if 'mykrobe' in config['analyses'] else mykrobe.OUTPUT_MYKROBE_REPORT_EMPTY),
        # Sequence typing
        report_mlst_warwick = sequence_typing.get_sequence_typing_report('mlst_warwick', config),
        report_mlst_pasteur = sequence_typing.get_sequence_typing_report('mlst_pasteur', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        report_citations = Path(config['working_dir'], core.OUTPUT_HTML_CITATIONS),
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name=config['sample_name'],
        output_dir=config['output_dir'],
        pipeline_info=config['pipeline'],
        input_dict=config['input'],
        input_type=config['input_type'],
        detection_method=config['detection_method'],
        citation_keys=config['citations']
    run:
        import datetime
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.components.pipelines.reportpipeline import ReportPipeline

        # Add the header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            sample_name=params.sample_name,
            date=datetime.datetime.now(),
            pipeline_version=params.pipeline_info['version'],
            input_files=ReportPipeline.format_input_string(params.input_dict),
            input_type=params.input_type,
            extra_data=[('Detection method', params.detection_method)],
            key_citation=params.citation_keys['main']
        ))

        # Add output sections
        report_structure = []
        ReportPipeline.add_content_scrubbing(
            report_structure, params.input_type, input.reports_scrubbing)
        ReportPipeline.add_content_trim_basic_qc(
            report_structure,params.input_type,input.reports_downsampling,input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        ReportPipeline.add_content_contamination_check(
            report_structure, params.input_type, input.reports_contamination, input.report_confindr)
        report_structure.extend([
            ('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]),
            ('Species identification', 'species', [Path(input.report_rmlst)]),
            ('<i>Shigella</i> serotyping', 'shigella_typing', [Path(x) for x in (
                input.report_shigeifinder, input.report_shigatyper, input.report_mykrobe)]),
            ('AMR detection', 'amr', [Path(x) for x in (
                input.report_amrfinder, input.report_resfinder4)]),
            ('Virulence characterization', 'viru', [Path(x) for x in (
                input.report_virulence, input.report_virulence_shiga)]),
            ('Plasmid replicon detection', 'plasmid', [Path(x) for x in (
                input.report_mob_suite, input.report_genomic_context)]),
            ('Sequence typing', 'st', [Path(x) for x in (
                input.report_mlst_warwick, input.report_mlst_pasteur, input.report_cgmlst)]),
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
        human_read_scrubbing.get_summaries(config),
        downsampling.get_summaries(config),
        trimming.get_summaries(config),
        Path(config['working_dir']) / quast.OUTPUT_QUAST_SUMMARY,
        contamination_check_kraken.get_summaries(config),
        confindr.get_summary(config),
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='rmlst') if 'rmlst' in config['analyses'] else [],
        # Shigella typing
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='serotype_o') if 'serotype_o' in config['analyses'] else [],
        # Shigella typing
        Path(config['working_dir']) / shigeifinder.OUTPUT_SHIGEIFINDER_SUMMARY if 'shigeifinder' in config['analyses'] else [],
        Path(config['working_dir']) / shigatyper.OUTPUT_SHIGATYPER_SUMMARY if 'shigatyper' in config['analyses'] else [],
        Path(config['working_dir']) / mykrobe.OUTPUT_MYKROBE_SUMMARY if 'mykrobe' in config['analyses'] else [],
         # Gene detection
        Path(config['working_dir']) / amrfinder.OUTPUT_AMRFINDER_SUMMARY if 'amrfinder' in config['analyses'] else [],
        Path(config['working_dir']) / resfinder4.OUTPUT_RESFINDER4_SUMMARY if 'resfinder4' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='virulencefinder_shiga') if 'virulencefinder' in config['analyses'] else [],
        # Plasmid characterization
        Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_SUMMARY if 'mob_suite' in config['analyses'] else [],
        # Sequence typing
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst_pasteur') if 'mlst_pasteur' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst_warwick') if 'mlst_warwick' in config['analyses'] else [],
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
        TSV_gd_coli = Path(config['working_dir']) / 'gene_detection' / 'virulencefinder' / 'metadata' / 'tsv.io' if 'virulencefinder' in config['analyses'] else [],
        INFORMS_gd_coli = Path(config['working_dir']) / 'gene_detection' / 'virulencefinder' / 'db_manager' / 'informs.io' if 'virulencefinder' in config['analyses'] else [],
        TSV_gd_shiga = Path(config['working_dir']) / 'gene_detection' / 'virulencefinder_shiga' / 'metadata' / 'tsv.io' if 'virulencefinder' in config['analyses'] else [],
        INFORMS_gd_shiga = Path(config['working_dir']) / 'gene_detection' / 'virulencefinder_shiga' / 'db_manager' / 'informs.io' if 'virulencefinder' in config['analyses'] else []
    output:
        TSV = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'input' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'input' / 'informs.io'
    run:
        mobsuite.collect_genomic_context_input(input, Path(output.TSV), Path(output.INFORMS))
