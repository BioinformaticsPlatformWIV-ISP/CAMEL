import shutil
from pathlib import Path

from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import trimming_illumina, trimming_ont, trimming, downsampling, \
    contamination_check_kraken, core, human_read_scrubbing
from camel.scripts.viralconsensuspipeline.snakefile import iterativemapping, refselection, preprocess, \
    multiallelicsites, nextclade3

#######################
# Included snakefiles #
#######################
include: core.SNAKEFILE_CORE
include: human_read_scrubbing.SNAKEFILE_SCRUBBING
include: downsampling.SNAKEFILE_DOWNSAMPLING
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: trimming_ont.SNAKEFILE_TRIMMING_ONT
include: refselection.SNAKEFILE_REF_SELECTION
include: preprocess.SNAKEFILE_PREPROCESS
include: iterativemapping.SNAKEFILE_ITERATIVE_MAPPING
include: nextclade3.SNAKEFILE_NEXTCLADE
include: multiallelicsites.SNAKEFILE_MULTI_ALLELIC


#########
# Rules #
#########
rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        HTML = config['output_report'],
        TSV = config['output_tabular'],
        TSV_stats = Path(config['working_dir']) / 'preprocess' / 'stats.tsv' if config['input_type'] != 'fasta' else []

rule link_fasta_to_iterative_mapping:
    """
    Selects the FASTA file used as a reference for read mapping.
    """
    input:
        FASTA = Path(config['working_dir']) / refselection.OUTPUT_REF_SELECTION_FASTA if config['fasta_ref'] is None else []
    output:
        FASTA = Path(config['working_dir']) / iterativemapping.INPUT_FASTA_REF
    params:
        input_type = config['input_type'],
        fasta_ref = config['fasta_ref']
    run:
        from camel.app.io.tooliofile import ToolIOFile
        if params.fasta_ref is not None:
            SnakemakeUtils.dump_object([ToolIOFile(Path(params.fasta_ref))], Path(output.FASTA))
        else:
            shutil.copyfile(input.FASTA, output.FASTA)

rule select_fasta_file:
    """
    This rule selects the fasta file to send to other workflows.
    """
    input:
        FASTA = Path(config['working_dir']) / iterativemapping.get_fasta(config)
    output:
        FASTA = Path(config['working_dir']) / 'fasta.io'
    shell:
        "cp {input.FASTA} {output.FASTA};"

rule report_create_command_section:
    """
    Creates the report section containing the tool commands.
    """
    input:
        INFORMS_scrubbing = human_read_scrubbing.get_command_informs(config),
        INFORMS_downsampling = downsampling.get_command_informs(config),
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        INFORMS_reference_selection = Path(config['working_dir']) / 'ref_selection' / 'mash_screen' / refselection.get_segments(
            Path(config['ref_selection']['db']))[0] / 'informs.io' if config['fasta_ref'] is None and config['input_type'] != 'fasta' and config['ref_selection'].get('db') is not None else [],
        INFORMS_preprocess = Path(config['working_dir']) / preprocess.OUTPUT_PRE_PROCESS_INFORMS if config['input_type'] != 'fasta' else [],
        INFORMS_iterative_mapping = Path(config['working_dir']) / iterativemapping.OUTPUT_ITERATIVE_MAPPING_INFORMS if config['input_type'] != 'fasta' else [],
        INFORMS_mash = Path(config['working_dir']) / 'nextclade' / 'subtype_determination' / 'mash' / 'informs.io' if config['nextclade'].get('db_mash') is not None and config['input_type'] != 'fasta' else [],
        INFORMS_nextclade = Path(config['working_dir']) / nextclade3.OUTPUT_NEXTCLADE_INFORMS if 'nextclade' in config['analyses'] else []
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        working_dir = config['working_dir']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.export_command_section(input, Path(output.HTML), params.working_dir)

rule report_combine_all:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        reports_scrubbing = human_read_scrubbing.get_reports(config),
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_reference_selection = Path(config['working_dir']) / (refselection.OUTPUT_REF_SELECTION_REPORT if config['fasta_ref'] is None and config['input_type'] != 'fasta' else refselection.OUTPUT_REF_SELECTION_REPORT_EMPTY),
        report_preprocess_ampligone = Path(config['working_dir']) / (preprocess.OUTPUT_PRE_PROCESS_AMPLIGONE_REPORT if 'ampligone' in config['analyses'] and config['input_type'] != 'fasta' else preprocess.OUTPUT_PRE_PROCESS_AMPLIGONE_REPORT_EMPTY),
        report_preprocess_clipping = Path(config['working_dir']) / (preprocess.OUTPUT_PRE_PROCESS_CLIPPING_REPORT if 'ampligone' in config['analyses'] and config['input_type'] != 'fasta' else preprocess.OUTPUT_PRE_PROCESS_CLIPPING_REPORT_EMPTY),
        report_preprocess = Path(config['working_dir']) / preprocess.OUTPUT_PRE_PROCESS_REPORT if config['input_type'] != 'fasta' else [],
        report_iterative_mapping = Path(config['working_dir']) / iterativemapping.OUTPUT_ITERATIVE_MAPPING_REPORT if config['input_type'] != 'fasta' else [],
        report_nexclade_subtype = Path(config['working_dir']) / (nextclade3.OUTPUT_NEXTCLADE_SUBTYPE_REPORT if (config['nextclade'].get('db') is None) and ('nextclade' in config['analyses']) and config['input_type'] != 'fasta' else nextclade3.OUTPUT_NEXTCLADE_SUBTYPE_REPORT_EMPTY),
        report_nextclade = Path(config['working_dir']) / (nextclade3.OUTPUT_NEXTCLADE_REPORT if 'nextclade' in config['analyses'] else nextclade3.OUTPUT_NEXTCLADE_REPORT_EMPTY),
        report_multi_allelic = Path(config['working_dir']) / multiallelicsites.OUTPUT_MULTI_ALLELIC_REPORT if config['input_type'] != 'fasta' else [],
        report_commands = Path(config['working_dir']) / 'report' / 'html-commands.io',
        report_citations = Path(config['working_dir'],core.OUTPUT_HTML_CITATIONS)
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        ref_genome = config['fasta_ref'],
        input_type = config['input_type'],
        config_input = config['input'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline']
    run:
        import datetime
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        # Add header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        extra_info = [('Reference genome', Path(params.ref_genome).name if params.ref_genome is not None else '-')]
        input_file_str = ', '.join(f['name'] for _, input_files in params.config_input.items() for f in input_files)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            params.sample_name,
            datetime.datetime.now(),
            params.pipeline_info['version'], input_file_str, input_type=params.input_type, extra_data=extra_info))

        # Other sections
        report_structure = []

        # Core sections (shared)
        ReportPipeline.add_content_scrubbing(
            report_structure, params.input_type, input.reports_scrubbing)
        ReportPipeline.add_content_trim_basic_qc(
            report_structure, params.input_type, input.reports_downsampling, input.reports_trimming)
        ReportPipeline.add_content_contamination_check(
            report_structure, params.input_type, input.reports_contamination, None)

        # Add output sections
        if params.input_type != 'fasta':
            report_structure.extend([
                ('Reference selection', 'ref_selection', [Path(input.report_reference_selection)]),
                ('Pre-processing', 'pre_process', [Path(x) for x in (
                    input.report_preprocess_ampligone, input.report_preprocess_clipping, input.report_preprocess)]),
                ('Consensus extraction', 'consensus', [Path(input.report_iterative_mapping)]),
                ('Multi-allelic sites', 'multi_allelic', [Path(input.report_multi_allelic)])
            ])
        report_structure.extend([
            ('Nextclade', 'nextclade', [Path(input.report_nexclade_subtype), Path(input.report_nextclade)]),
            ('Commands', 'commands', [Path(input.report_commands)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
        ])
        SnakePipelineUtils.add_report_content(report, report_structure)

rule summary_combine_all:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        Path(config['working_dir'], core.OUTPUT_TSV_SUMMARY_INIT),
        human_read_scrubbing.get_summaries(config),
        downsampling.get_summaries(config),
        trimming.get_summaries(config),
        contamination_check_kraken.get_summaries(config),
        Path(config['working_dir']) / refselection.OUTPUT_REF_SELECTION_SUMMARY if config['fasta_ref'] is None and config['input_type'] != 'fasta' else [],
        Path(config['working_dir']) / preprocess.OUTPUT_PRE_PROCESS_SUMMARY if config['input_type'] != 'fasta' else [],
        Path(config['working_dir']) / iterativemapping.OUTPUT_ITERATIVE_MAPPING_SUMMARY if config['input_type'] != 'fasta' else [],
        Path(config['working_dir']) / multiallelicsites.OUTPUT_MULTI_ALLELIC_SUMMARY if config['input_type'] != 'fasta' else [],
        Path(config['working_dir']) / nextclade3.OUTPUT_NEXTCLADE_SUMMARY if 'nextclade' in config['analyses'] else []
    output:
        config.get('output_tabular')
    run:
        with open(output[0], 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
