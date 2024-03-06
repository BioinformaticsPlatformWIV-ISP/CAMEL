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
        OUT = Path(config['working_dir']) / human_read_scrubbing.OUTPUT_SCRUBBING_FASTA if 'fasta' in config['input'] else Path(config['working_dir']) / human_read_scrubbing.OUTPUT_SCRUBBING_FASTQ


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
        INFORMS_scrubbing = Path(config['working_dir']) / human_read_scrubbing.OUTPUT_SCRUBBING_INFORMS
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
        report_scrubbing =  Path(config['working_dir']) / human_read_scrubbing.OUTPUT_SCRUBBING_REPORT,
        # Report
        # report_citations = rules.report_pickle_citations.output.HTML,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        config_input = config['input'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        citation_keys = config['citations'],
        input_type = config['input_type']
    run:
        import datetime
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        # Add header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            params.sample_name,
            datetime.datetime.now(),
            params.pipeline_info['version'],
            ', '.join(
                input_file['name'] for _, input_files in params.config_input.items() for input_file in input_files),
            params.input_type,
            params.citation_keys['main']
        ))

        # Add report content
        report_structure = [
            ('Human read scrubbing', 'scrub', [Path(input.report_scrubbing)]),
            # ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ]
        SnakePipelineUtils.add_report_content(report, report_structure)
