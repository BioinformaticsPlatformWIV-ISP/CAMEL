from pathlib import Path

from camel.resources.snakefile import trimming_illumina, downsampling, trimming_ont, trimming

#######################
# Included snakefiles #
#######################
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: trimming_ont.SNAKEFILE_TRIMMING_ONT
include: downsampling.SNAKEFILE_DOWNSAMPLING

#########
# Rules #
#########
rule all:
    """
    Ensures that the required output files are generated.
    """
    input:
        config['output_report'],
        config['output_summary']

rule report_create:
    """
    Creates the output HTML report.
    """
    input:
        report_downsampling = Path(config['working_dir']) /downsampling.OUTPUT_DOWNSAMPLING_REPORT,
        report_trimming = trimming.get_trimming_report(config)
    output:
        HTML = config['output_report']
    params:
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline']
    run:
        import datetime
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        # Add the header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            'test_sample',
            datetime.datetime.now(),
            params.pipeline_info['version'], 'todo', []))

        # Add output sections
        report_structure = [
            ('Read trimming and basic QC', 'trim', [Path(input.report_downsampling), Path(input.report_trimming)]),
        ]
        SnakePipelineUtils.add_report_content(report, report_structure)
