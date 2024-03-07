import shutil
from pathlib import Path

from camel.resources.snakefile import trimming_illumina, downsampling, trimming_ont, trimming, quast, \
    contamination_check_kraken, quality_checks, confindr, gene_detection, assembly, core, assembly_flye, \
    polish_assembly_long

#######################
# Included snakefiles #
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

#########s
# Rules #
#########
rule all:
    """
    Ensures that the required output files are generated.
    """
    input:
        config['output_report'],
        config['output_tabular']

rule copy_assemblies_to_medaka_input:
    input:
        FASTA_flye = Path(config['working_dir']) / assembly_flye.OUTPUT_ASSEMBLY_FASTA
    output:
        FASTA_medaka_flye = str(Path(config['working_dir']) / polish_assembly_long.INPUT_ASSEMBLY_FASTA).format(assembly_type='flye')
    run:
        shutil.copyfile(input.FASTA_flye, output.FASTA_medaka_flye)

rule report_create_command_section:
    """
    Creates the report section containing the tool commands.
    """
    input:
        INFORMS_downsampling = downsampling.get_command_informs(config),
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_assembly = assembly.get_command_informs(config),
        INFORMS_quast = Path(config['working_dir'], quast.OUTPUT_QUAST_INFORMS),
        INFORMS_busco = Path(config['working_dir'], quast.OUTPUT_BUSCO_INFORMS),
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        INFORMS_confindr = confindr.get_command_informs(config),
        INFORMS_ncbi_amr = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else []
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.export_command_section(input, Path(output.HTML), params.dir_)

rule report_create:
    """
    Creates the output HTML report.
    """
    input:
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        report_quast = Path(config['working_dir'], quast.OUTPUT_QUAST_REPORT),
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_confindr = confindr.get_report(config),
        report_adv_qc = Path(config['working_dir']) / str(quality_checks.OUTPUT_QUALITY_CHECKS_REPORT).format(
            input_type=config['input_type']),
        report_ncbi_amr = gene_detection.get_gene_detection_report('ncbi_amr', config),
        report_commands= rules.report_create_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        output_dir = config['output_dir'],
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
        ReportPipeline.add_content_trim_basic_qc(
            report_structure, params.input_type, input.reports_downsampling, input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        ReportPipeline.add_content_contamination_check(
            report_structure, params.input_type, input.reports_contamination, input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))
        report_structure.append(('Gene detection', 'gene_detection', [Path(input.report_ncbi_amr)]))
        report_structure.append(('Commands', 'commands', [Path(input.report_commands)]))

        # Export the report
        SnakePipelineUtils.add_report_content(report, report_structure)

rule summary_combine:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        Path(config['working_dir'], core.OUTPUT_TSV_SUMMARY_INIT),
        downsampling.get_summaries(config),
        trimming.get_summaries(config),
        Path(config['working_dir'], quast.OUTPUT_QUAST_SUMMARY),
        contamination_check_kraken.get_summaries(config),
        confindr.get_summary(config),
        Path(config['working_dir'], quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY)
    output:
        TSV = config.get('output_tabular')
    run:
        with open(output.TSV, 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
