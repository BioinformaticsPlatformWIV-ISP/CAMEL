from pathlib import Path

from camel.resources.snakefile import human_read_scrubbing

#######################
# Included Snakefiles #
#######################
include: human_read_scrubbing.SNAKEFILE_SCRUBBING

#########
# Rules #
#########
rule all:
    input:
        HTML = config['output_report'],
        OUT = Path(config['working_dir']) / 'human_read_scrubbing' / 'fasta' / 'output' / 'fasta.io' if 'fasta' in config['input'] else Path(config['working_dir']) / 'human_read_scrubbing' / ('fastq_pe' if 'fastq_pe' in config['input'] else 'fastq_se') / 'output' / 'fastq.io'


rule link_scrubbing_input:
    """
    Creates the FASTQ input for the human read scrubbing step.
    """
    output:
        FASTQ = Path(config['working_dir']) / human_read_scrubbing.INPUT_SCRUBBING_FASTQ if not 'fasta' in config['input'] else [],
        FASTA = Path(config['working_dir']) / human_read_scrubbing.INPUT_SCRUBBING_FASTA if 'fasta' in config['input'] else []
    run:
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        from camel.app.io.tooliofile import ToolIOFile
        if 'fasta' in config['input']:
            SnakemakeUtils.dump_object([ToolIOFile(Path(config['input']['fasta'][0]['path']))], Path(output.FASTA))
        elif 'fastq_pe' in config['input']:
            SnakemakeUtils.dump_object([ToolIOFile(Path(f['path'])) for f in config['input']['fastq_pe']], Path(output.FASTQ))
        elif 'fastq_se' in config['input']:
            SnakemakeUtils.dump_object([ToolIOFile(Path(f['path'])) for f in config['input']['fastq_se']], Path(output.FASTQ))


rule report_command_section:
    """
    Creates the report section containing the tool commands.
    """
    input:
        INFORMS_scrubbing = human_read_scrubbing.get_command_informs(config)
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
        # Report
        # report_citations = rules.report_pickle_citations.output.HTML,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        input_dict = config['input'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        citation_keys = config['citations'],
        input_type = config['input_type']
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
        ))

        # Add report content
        report_structure = []
        ReportPipeline.add_content_scrubbing(
            report_structure, params.input_type, input.reports_scrubbing)
        report_structure.extend([('Commands', 'commands', [Path(input.report_commands)])])
        SnakePipelineUtils.add_report_content(report, report_structure)
