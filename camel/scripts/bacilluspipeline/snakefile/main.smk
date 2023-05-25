import shutil
from pathlib import Path

from camel.resources.snakefile import trimming, trimming_illumina, assembly_spades, assembly_canu, \
    quality_checks, contamination_check_kraken, sequence_typing, downsampling, amrfinder, trimming_ont, gene_detection
from camel.scripts.bacilluspipeline.snakefile import btyper

#######################
# Included Snakefiles #
#######################
include: downsampling.SNAKEFILE_DOWNSAMPLING
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: trimming_ont.SNAKEFILE_TRIMMING_ONT
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: assembly_spades.SNAKEFILE_ASSEMBLY_SPADES
include: assembly_canu.SNAKEFILE_ASSEMBLY_CANU
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING
include: btyper.SNAKEFILE_BTYPER
include: amrfinder.SNAKEFILE_AMRFINDER
include: gene_detection.SNAKEFILE_GENE_DETECTION

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
        IO_FASTQ= Path(config['working_dir']) / 'fq_dict.io'

rule link_downsampling_input:
    """
    Creates the FASTQ input for the downsampling step. 
    """
    output:
        FASTQ = Path(config['working_dir']) / downsampling.INPUT_DOWNSAMPLING_FASTQ
    params:
        config_input=config['input'],
        read_type=config['read_type']
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils

        if params.read_type == 'illumina':
            SnakemakeUtils.dump_object([ToolIOFile(Path(x['path'])) for x in config['input']['fastq_pe']],
                Path(output.FASTQ))
        elif params.read_type == 'nanopore':
            SnakemakeUtils.dump_object([ToolIOFile(Path(x['path'])) for x in config['input']['fastq_se']],
                Path(output.FASTQ))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

rule link_downsampling_to_trimming_workflows:
    """
    Links the downsampling output to the input of the ONT trimming workflow.  
    """
    input:
        FASTQ = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_FASTQ
    output:
        FASTQ_ilmn = Path(config['working_dir']) / trimming_illumina.INPUT_TRIMMOMATIC_FASTQ if config['read_type'] == 'illumina' else [],
        FASTQ_ont = Path(config['working_dir']) / trimming_ont.INPUT_ONT_FASTQ if config['read_type'] == 'nanopore' else []
    params:
        read_type=config['read_type']
    run:
        if params.read_type == 'nanopore':
            shutil.copyfile(Path(input.FASTQ),Path(output.FASTQ_ont))
        elif params.read_type == 'illumina':
            shutil.copyfile(Path(input.FASTQ),Path(output.FASTQ_ilmn))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

rule select_fastq_to_io:
    """
    Creates an IO object with the trimmed FASTQ files.
    Other workflows such as Kraken or de novo assembly rely on this dictionary to get input files (PE or SE).
    """
    input:
        FASTQ_ilmn = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_DICT if config['read_type'] == 'illumina' else [],
        FASTQ_ont = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_DICT if config['read_type'] == 'nanopore' else []
    output:
        IO_FASTQ = Path(config['working_dir']) / 'fq_dict.io'
    params:
        read_type = config['read_type']
    run:
        if params.read_type == 'illumina':
            shutil.copyfile(Path(input.FASTQ_ilmn), Path(output.IO_FASTQ))
        elif params.read_type == 'nanopore':
            shutil.copyfile(Path(input.FASTQ_ont), Path(output.IO_FASTQ))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

# For some reason, I cannot include the lines in the rule select_fasta_to_tools
rule select_fasta_to_genedetection:
    """
        This rules links the output of the assembly workflow to the other workflows.
    """
    input:
        FASTA_spades = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA if config['read_type'] == 'illumina' else [],
        FASTA_canu = Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_FASTA if config['read_type'] == 'nanopore' else []
    output:
        FASTA_genedetection = Path(config['working_dir']) / gene_detection.INPUT_GENE_DETECTION_FASTA
    params:
        read_type=config['read_type']
    run:
        if params.read_type == 'nanopore':
            shutil.copyfile(Path(input.FASTA_canu), Path(output.FASTA_genedetection))
        elif params.read_type == 'illumina':
            shutil.copyfile(Path(input.FASTA_spades),Path(output.FASTA_genedetection))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

rule select_fasta_to_tools:
    """
    This rules links the output of the assembly workflow to the other workflows.
    """
    input:
        FASTA_spades = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA if config['read_type'] == 'illumina' else [],
        FASTA_canu = Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_FASTA if config['read_type'] == 'nanopore' else []
    output:
        FASTA_btyper = Path(config['working_dir']) / btyper.INPUT_BTYPER_FASTA,
        FASTA_amrfinder = Path(config['working_dir']) / amrfinder.INPUT_RESFINDER_FASTA,
    params:
        read_type=config['read_type']
    run:
        if params.read_type == 'nanopore':
            shutil.copyfile(Path(input.FASTA_canu), Path(output.FASTA_btyper))
            shutil.copyfile(Path(input.FASTA_canu), Path(output.FASTA_amrfinder))
        elif params.read_type == 'illumina':
            shutil.copyfile(Path(input.FASTA_spades),Path(output.FASTA_btyper))
            shutil.copyfile(Path(input.FASTA_spades),Path(output.FASTA_amrfinder))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

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
        INFORMS_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_INFORMS,
        INFORMS_trimming = trimming.get_trimming_command_informs(config),
        INFORMS_assembly = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_INFORMS if config['read_type'] == 'illumina' else Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_INFORMS,
        INFORMS_assembly_filt = Path(config['working_dir']) / 'assembly_spades' / 'filtering' / 'informs.io' if config['read_type'] == 'illumina' else [],
        INFORMS_kraken = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS if 'kraken' in config['analyses'] else [],
        INFORMS_gmo = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='gmo') if 'gmo' in config['analyses'] else[],
        INFORMS_plasmid = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='plasmid') if 'plasmid' in config['analyses'] else[],
        INFORMS_mapping = quality_checks.get_mapping_rate_informs(config) if config['read_type'] == 'illumina' else [],
        INFORMS_depth = quality_checks.get_depth_informs(config) if config['read_type'] == 'illumina' else [],
        INFORMS_btyper = Path(config['working_dir']) / str(btyper.OUTPUT_INFORMS_BTYPER).format(scheme='btyper_typing') if 'btyper' in config['analyses'] else [],
        INFORMS_amrfinder = Path(config['working_dir']) / str(amrfinder.OUTPUT_AMRFINDER_INFORMS) if 'amrfinder' in config['analyses'] else []
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
        report_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_REPORT,
        report_trimming = trimming.get_trimming_report(config),
        report_assembly = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_REPORT if config['read_type'] == 'illumina' else Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_REPORT,
        report_kraken = Path(config['working_dir']) / (contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT if 'kraken' in config['analyses'] else contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY),
        report_adv_qc = Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT if config['read_type'] == 'illumina' else [],
        report_gmo = gene_detection.get_gene_detection_report('gmo',config),
        report_plasmid = gene_detection.get_gene_detection_report('plasmid',config),
        report_mlst=sequence_typing.get_sequence_typing_report('mlst',config),
        report_cgmlst=sequence_typing.get_sequence_typing_report('cgmlst',config),
        report_btyper = Path(config['working_dir']) / btyper.OUTPUT_BTYPER_REPORT,
        report_amrfinder = Path(config['working_dir']) / amrfinder.OUTPUT_AMRFINDER_REPORT,
        report_citations = rules.report_pickle_citations.output.HTML,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        config_input = config['input'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        detection_method = config['detection_method'],
        citation_keys = config['citations']
    run:
        import datetime
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        # Add header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            params.sample_name,
            datetime.datetime.now(),
            params.pipeline_info['version'], ', '.join(
                input_file['name'] for _, input_files in params.config_input.items() for input_file in input_files),
            [('Detection method', params.detection_method)],
            params.citation_keys['main']
        ))

        # Add content
        report_structure = [
            ('Read trimming and basic QC', 'trim', [Path(input.report_downsampling), Path(input.report_trimming)]),
            ('Assembly', 'assem', [Path(input.report_assembly)]),
            ('Advanced QC', 'adv_qc', [Path(x) for x in (input.report_kraken, input.report_adv_qc)] if config['read_type'] == 'illumina' else [Path(input.report_kraken)]),
            ('GMO detection', 'gmo', [Path(input.report_gmo)]),
            ('Plasmid detection', 'plasmid', [Path(input.report_plasmid)]),
            ('Sequence typing', 'st', [Path(x) for x in (input.report_mlst, input.report_cgmlst)]),
            ('BTyper results', 'btyper', [Path(input.report_btyper)]),
            ('AMRFinder results', 'amrfinder', [Path(input.report_amrfinder)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ]
        SnakePipelineUtils.add_report_content(report, report_structure)

rule summary_init:
    """
    Initializes the summary output file.
    """
    output:
        TSV = Path(config['working_dir']) / 'summary' / 'summary-init.tsv'
    run:
        import datetime
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        analysis_date = datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT)
        input_filenames = ', '.join(
            input_file['name'] for _, input_files in config['input'].items() for input_file in input_files)
        with open(output.TSV, 'w') as handle:
            for kv_pair in [
                ('pipeline_name', config['pipeline']['name']),
                ('pipeline_version', config['pipeline']['version']),
                ('sample', config['sample_name']),
                ('input_files', input_filenames),
                ('analysis_date', analysis_date),
                ('detection_method', config['detection_method'])]:
                handle.write('\t'.join(kv_pair))
                handle.write('\n')

rule summary_combine_all:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        rules.summary_init.output.TSV,
        Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_SUMMARY,
        trimming.get_trimming_summary(config),
        Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_SUMMARY if config['read_type'] == 'illumina' else Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_SUMMARY,
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY if config['read_type'] == 'illumina' else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='gmo') if 'gmo' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='plasmid') if 'plasmid' in config['analyses'] else [],
        Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_SUMMARY if 'kraken' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        Path(config['working_dir']) / btyper.OUTPUT_BTYPER_SUMMARY if 'btyper' in config['analyses'] else [],
        Path(config['working_dir']) / amrfinder.OUTPUT_AMRFINDER_SUMMARY if 'amrfinder' in config['analyses'] else []
    output:
        config.get('output_tabular')
    run:
        with open(output[0], 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
