from pathlib import Path

from camel.resources.snakefile import trimming_illumina, gene_detection, trimming, contamination_check_kraken, \
    quality_checks, sequence_typing, lrefinder, downsampling, confindr, quast, core, assembly, amrfinder, resfinder4, \
    bacmet, mobsuite

#######################
# Included Snakefiles #
#######################
include: core.SNAKEFILE_CORE
include: downsampling.SNAKEFILE_DOWNSAMPLING
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: confindr.SNAKEFILE_CONFINDR
include: assembly.SNAKEFILE_ASSEMBLY
include: quast.SNAKEFILE_QUAST
include: lrefinder.SNAKEFILE_LREFINDER
include: amrfinder.SNAKEFILE_AMRFINDER
include: resfinder4.SNAKEFILE_RESFINDER4
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: mobsuite.SNAKEFILE_MOB_SUITE
include: bacmet.SNAKEFILE_BACMET
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING


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
        INFORMS_downsampling = downsampling.get_command_informs(config),
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_assembly = assembly.get_command_informs(config),
        INFORMS_quast = Path(config['working_dir']) / quast.OUTPUT_QUAST_INFORMS,
        INFORMS_busco = Path(config['working_dir']) / quast.OUTPUT_BUSCO_INFORMS,
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        INFORMS_confindr = Path(config['working_dir']) / confindr.OUTPUT_CONFINDR_INFORMS if 'confindr' in config['analyses'] else [],
        INFORMS_assembly_map = assembly.get_qc_informs(config, config['input_type']),
        INFORMS_lrefinder= Path(config['working_dir']) / lrefinder.OUTPUT_LREFINDER_INFORMS if 'lrefinder' in config['analyses'] else [],
        INFORMS_amrfinder = Path(config['working_dir']) / amrfinder.OUTPUT_AMRFINDER_INFORMS if 'amrfinder' in config['analyses'] else [],
        INFORMS_resfinder4 = Path(config['working_dir']) / resfinder4.OUTPUT_RESFINDER4_INFORMS if 'resfinder4' in config['analyses'] else [],
        INFORMS_mob_suite = Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_INFORMS if 'mob_suite' in config['analyses'] else [],
        INFORMS_virulencefinder = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        INFORMS_vfdb_core = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        IFNORMS_bacmet = Path(config['working_dir']) / bacmet.OUTPUT_BACMET_INFORMS if 'bacmet' in config['analyses'] else [],
        INFORMS_prodigal = Path(config['working_dir']) / bacmet.OUTPUT_PRODIGAL_INFORMS if 'bacmet' in config['analyses'] else [],
        INFORMS_mlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='mlst') if 'mlst' in config['analyses'] else []
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        working_dir = config['working_dir']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        informs = []
        for content in [SnakemakeUtils.load_object(Path(io)) for io in input]:
            if type(content) is dict:
                informs.append(content)
            elif type(content) is list:
                informs.extend(content)
        section = SnakePipelineUtils.create_commands_section(informs, params.working_dir)
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule report_combine_all:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        report_quast = Path(config['working_dir']) /quast.OUTPUT_QUAST_REPORT,
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_confindr = Path(config['working_dir']) / (confindr.OUTPUT_CONFINDR_REPORT if 'confindr' in config['analyses'] else confindr.OUTPUT_CONFINDR_REPORT_EMPTY),
        report_adv_qc = Path(config['working_dir']) / str(quality_checks.OUTPUT_QUALITY_CHECKS_REPORT).format(input_type=config['input_type']),
        # AMR detection
        report_lrefinder = Path(config['working_dir']) / (lrefinder.OUTPUT_LREFINDER_REPORT if 'lrefinder' in config['analyses'] else lrefinder.OUTPUT_LREFINDER_REPORT_EMPTY),
        report_amrfinder = Path(config['working_dir']) / (amrfinder.OUTPUT_AMRFINDER_REPORT if config['analyses'] else amrfinder.OUTPUT_AMRFINDER_REPORT_EMPTY),
        report_resfinder4 = Path(config['working_dir']) / (resfinder4.OUTPUT_RESFINDER4_REPORT if config['analyses'] else resfinder4.OUTPUT_RESFINDER4_REPORT_EMPTY),
        # Virulence gene detection
        report_virulencefinder = gene_detection.get_gene_detection_report('virulencefinder', config),
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        # Plasmid characterization
        report_plasmidfinder = gene_detection.get_gene_detection_report('plasmidfinder', config),
        report_mob_suite = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_REPORT_EMPTY),
        report_genomic_context = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'html.io',
        # BacMet
        report_prodigal=Path(config['working_dir']) / (bacmet.OUTPUT_PRODIGAL_REPORT if 'bacmet' in config['analyses'] else bacmet.OUTPUT_PRODIGAL_REPORT_EMPTY),
        report_bacmet=Path(config['working_dir']) / (bacmet.OUTPUT_BACMET_REPORT if 'bacmet' in config['analyses'] else bacmet.OUTPUT_BACMET_REPORT_EMPTY),
        # Typing
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst',config),
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config),
        report_mlst_bezdicek = sequence_typing.get_sequence_typing_report('mlst_bezdicek', config),
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
        species = config['selected_species'],
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
        report.add_html_object(SnakePipelineUtils.create_input_section(
            sample_name=params.sample_name,
            date=datetime.datetime.now(),
            pipeline_version=params.pipeline_info['version'],
            input_files=ReportPipeline.format_input_string(params.input_dict),
            input_type=params.input_type,
            extra_data=[('Selected species', f'<i>{params.species}</i>')],
            key_citation=params.citation_keys['main']
        ))

        # Add report content
        report_structure = []
        ReportPipeline.add_content_trim_basic_qc(
            report_structure,params.input_type,input.reports_downsampling,input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        ReportPipeline.add_content_contamination_check(
            report_structure,params.input_type,input.reports_contamination,input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))

        # Typing (additional MLST scheme for E. faecium)
        if params.species == 'faecalis':
            reports_typing = (input.report_rmlst, input.report_mlst, input.report_cgmlst)
        else:
            reports_typing = (input.report_rmlst, input.report_mlst, input.report_mlst_bezdicek, input.report_cgmlst)

        report_structure.extend([
            ('AMR detection', 'amr', [Path(x) for x in (
                input.report_lrefinder, input.report_amrfinder, input.report_resfinder4)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (
                input.report_virulencefinder, input.report_vfdb_core)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_plasmidfinder, input.report_mob_suite, input.report_genomic_context)]),
            ('Biocide and metal resistance', 'bacmet', [Path(input.report_prodigal), Path(input.report_bacmet)]),
            ('Sequence typing', 'st', [Path(x) for x in reports_typing]),
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
        downsampling.get_summaries(config),
        trimming.get_summaries(config),
        Path(config['working_dir']) / quast.OUTPUT_QUAST_SUMMARY,
        contamination_check_kraken.get_summaries(config),
        confindr.get_summary(config),
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        # AMR detection
        Path(config['working_dir']) / lrefinder.OUTPUT_LREFINDER_SUMMARY if 'lrefinder' in config['analyses'] else [],
        Path(config['working_dir']) / amrfinder.OUTPUT_AMRFINDER_SUMMARY if 'amrfinder' in config['analyses'] else [],
        Path(config['working_dir']) / resfinder4.OUTPUT_RESFINDER4_SUMMARY if 'resfinder4' in config['analyses'] else [],
        # Virulence detection
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        # Plasmid characterization
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else [],
        Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_SUMMARY if 'mob_suite' in config['analyses'] else [],
        # BacMet
        Path(config['working_dir']) / bacmet.OUTPUT_BACMET_SUMMARY if 'bacmet' in config['analyses'] else [],
        # Sequence typing
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='rmlst') if 'rmlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else []
    output:
        TSV = config['output_tabular']
    run:
        with open(output.TSV, 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
