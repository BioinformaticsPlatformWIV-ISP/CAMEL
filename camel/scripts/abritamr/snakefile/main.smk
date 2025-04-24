from pathlib import Path
from camel.resources.snakefile import abritamr, assembly, core

#######################
# Included Snakefiles #
#######################
include: core.SNAKEFILE_CORE
include: abritamr.SNAKEFILE_ABRITAMR
include: assembly.SNAKEFILE_ASSEMBLY

#########
# Rules #
#########
rule all:
    input:
        HTML = config['output_report'],
        TSV = config['output_tabular']


rule report_pickle_citations:
    """
    This rule creates a pickle with a report section containing the citations.
    """
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-citations.io'
    params:
        citation_keys = config['citations']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        section = SnakePipelineUtils.create_citations_section(
            params.citation_keys['other'], params.citation_keys['main'])
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule report_command_section:
    input:
        INFORMS_abritamr_run =  Path(config['working_dir']) / str(abritamr.OUTPUT_ABRITAMR_RUN_INFORMS),
        INFORMS_abritamr_report=  Path(config['working_dir']) / str(abritamr.OUTPUT_REPORT_ABRITAMR_INFORMS)
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
        report_abritamr = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_REPORT,
        # Report
        report_citations = rules.report_pickle_citations.output.HTML,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        input_dict = config['input'],
        input_type = config['input_type'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        citation_keys = config['citations']
    run:
        import datetime
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        # Add header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            sample_name=params.sample_name,
            date=datetime.datetime.now(),
            pipeline_version=params.pipeline_info['version'],
            input_files=ReportPipeline.format_input_string(params.input_dict),
            input_type=params.input_type,
            key_citation=params.citation_keys['main'],
            warning=False
        ))

        # Add report content
        report_structure = [
            ('AbriTAMR', 'abrit', [Path(input.report_abritamr)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
                            ]
        SnakePipelineUtils.add_report_content(report, report_structure)

# rule summary_init:
#     """
#     Initializes the summary output file.
#     """
#     output:
#         TSV = Path(config['working_dir']) / 'summary' / 'summary-init.tsv',
#     run:
#         import datetime
#         from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
#         analysis_date = datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT)
#         input_filenames = ', '.join(
#             input_file['name'] for _, input_files in config['input'].items() for input_file in input_files)
#         with open(output.TSV, 'w') as handle:
#             for kv_pair in [
#                 ('pipeline_name', config['pipeline']['name']),
#                 ('pipeline_version', config['pipeline']['version']),
#                 ('sample', config['sample_name']),
#                 ('input_files', input_filenames),
#                 ('analysis_date', analysis_date)]:
#                 handle.write('\t'.join(kv_pair))
#                 handle.write('\n')
#                 json_dict[kv_pair[0]] = kv_pair[1]
#
# rule summary_combine_all:
#     """
#     In this rule all summary files are combined into a complete summary output file.
#     """
#     input:
#         rules.summary_init.output.TSV,
#         Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_SUMMARY,
#     output:
#         config.get('output_tabular')
#     run:
#         with open(output[0], 'w') as handle_out:
#             for summary_input in input:
#                 with open(summary_input) as handle_in:
#                     handle_out.write(handle_in.read())


rule summary_combine_all:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_SUMMARY
    output:
        config.get('output_tabular')
    run:
        with open(output[0],'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
