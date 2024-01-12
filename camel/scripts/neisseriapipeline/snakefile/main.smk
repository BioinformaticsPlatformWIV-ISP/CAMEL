from pathlib import Path

from camel.app.camel import Camel
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import trimming, trimming_illumina, quality_checks, \
    contamination_check_kraken, gene_detection, sequence_typing, downsampling, confindr, quast, core, trimming_ont, \
    assembly
from camel.scripts.neisseriapipeline.snakefile import serogroup_determination, gmats, mendevar


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
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING
include: gmats.SNAKEFILE_GMATS
include: mendevar.SNAKEFILE_MENDEVAR
include: serogroup_determination.SNAKEFILE_SEROGROUP_DETERMINATION

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
        INFORMS_assembly_map = assembly.get_qc_informs(config, config['input_type']),
        INFORMS_resfinder = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='resfinder') if 'resfinder' in config['analyses'] else [],
        INFORMS_ncbi_amr = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
        INFORMS_rmlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='rmlst') if 'rmlst' in config['analyses'] else [],
        INFORMS_mlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        INFORMS_rplf = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='rplf') if 'rplf' in config['analyses'] else [],
        INFORMS_bast = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='bast') if 'bast' in config['analyses'] else [],
        INFORMS_pora = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='pora') if 'pora' in config['analyses'] else [],
        INFORMS_porb = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='porb') if 'porb' in config['analyses'] else [],
        INFORMS_feta = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='feta') if 'feta' in config['analyses'] else [],
        INFORMS_amr = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='resistance_genes') if 'resistance_genes' in config['analyses'] else [],
        INFORMS_vaccine = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='vaccine_targets') if 'vaccine_targets' in config['analyses'] else [],
        INFORMS_fhbp = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='fhbp') if 'fhbp' in config['analyses'] else [],
        INFORMS_cgmlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        INFORMS_serogroup = Path(config['working_dir']) / serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_INFORMS if 'serogroup' in config['analyses'] else []
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.export_command_section(input, Path(output.HTML), Path(params.dir_))

rule neisseria_additional_resistance_gene_metadata:
    """
    This rule is used to add resistance gene metadata for penA and rpoB genes.
    The data is parsed from the PubMLST webpage.
    """
    input:
        hits = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_HITS).format(scheme='resistance_genes', locus_type='DNA', detection_method=config['detection_method']),
        VAL_HTML = sequence_typing.get_sequence_typing_report('resistance_genes', config),
        INFORMS_scheme = Path(config['working_dir']) / 'typing' / 'resistance_genes' / 'informs-locus_set.io'
    output:
        VAL_HTML = Path(config['working_dir']) / 'typing' / 'resistance_genes' / 'metadata' / 'report.html'
    params:
        working_dir = Path(config['working_dir']) / 'typing' / 'resistance_genes' / 'metadata',
        loci='penA, rpoB'
    run:
        from camel.app.pipeline.step import Step
        from camel.app.tools.pipelines.neisseria.resistancemetadataextractor import ResistanceMetadataExtractor
        extractor = ResistanceMetadataExtractor(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(extractor, input)
        step = Step(str(rule), extractor, Camel.get_instance(), params.working_dir, config)
        extractor.update_parameters(loci=params.loci)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(extractor, output)

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
        report_adv_qc=Path(config['working_dir']) / str(quality_checks.OUTPUT_QUALITY_CHECKS_REPORT).format(
            input_type=config['input_type']),
        report_resfinder = gene_detection.get_gene_detection_report('resfinder', config),
        report_ncbi_amr = gene_detection.get_gene_detection_report('ncbi_amr', config),
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        report_pora = sequence_typing.get_sequence_typing_report('pora', config),
        report_porb = sequence_typing.get_sequence_typing_report('porb', config),
        report_feta = sequence_typing.get_sequence_typing_report('feta', config),
        report_rplf = sequence_typing.get_sequence_typing_report('rplf', config),
        report_vaccine_targets = sequence_typing.get_sequence_typing_report('vaccine_targets', config),
        report_resistance_genes = rules.neisseria_additional_resistance_gene_metadata.output.VAL_HTML if 'resistance_genes' in config['analyses'] else Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_REPORT_EMPTY).format(scheme='resistance_genes'),
        report_fhbp = sequence_typing.get_sequence_typing_report('fhbp', config),
        report_bast = sequence_typing.get_sequence_typing_report('bast', config),
        report_gmats = Path(config['working_dir']) / (gmats.OUTPUT_GMATS_REPORT if 'gmats' in config['analyses'] else gmats.OUTPUT_GMATS_REPORT_EMPTY),
        report_mendevar = Path(config['working_dir']) / (mendevar.OUTPUT_MENDEVAR_REPORT if 'mendevar' in config['analyses'] else mendevar.OUTPUT_MENDEVAR_REPORT_EMPTY),
        report_serogroup = Path(config['working_dir']) / (serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_REPORT if 'serogroup' in config['analyses'] else serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_REPORT_EMPTY),
        report_serogroup_legacy = Path(config['working_dir']) / (serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_LEGACY_REPORT if 'serogroup' in config['analyses'] else serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_LEGACY_REPORT_EMPTY),
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
            key_citation=params.citation_keys['main']
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
            ('Resistance characterization', 'res', [Path(x) for x in (input.report_resfinder, input.report_ncbi_amr)]),
            ('Sequence typing', 'st', [Path(x) for x in (
                input.report_rmlst, input.report_mlst, input.report_rplf, input.report_pora, input.report_porb,
                input.report_feta, input.report_resistance_genes, input.report_vaccine_targets, input.report_fhbp,
                input.report_cgmlst)]),
            ('Antigen typing', 'at', [Path(x) for x in (input.report_bast, input.report_gmats, input.report_mendevar)]),
            ('Serogroup determination', 'serogroup', [Path(
                input.report_serogroup), Path(input.report_serogroup_legacy)]),
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
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='resfinder') if 'resfinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='rmlst') if 'rmlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='rplf') if 'rplf' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='bast') if 'bast' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='pora') if 'pora' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='porb') if 'porb' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='feta') if 'feta' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='fhbp') if 'fhbp' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='resistance_genes') if 'resistance_genes' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='vaccine_targets') if 'vaccine_targets' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        Path(config['working_dir']) / gmats.OUTPUT_GMATS_SUMMARY if 'gmats' in config['analyses'] else [],
        Path(config['working_dir']) / serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_SUMMARY if 'serogroup' in config['analyses'] else []
    output:
        TSV = config.get('output_tabular')
    run:
        with open(output.TSV, 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
