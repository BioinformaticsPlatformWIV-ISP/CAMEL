from pathlib import Path

from camel.resources.snakefile import trimming_illumina, assembly_spades, gene_detection, trimming, \
    contamination_check_kraken, quality_checks, sequence_typing
from camel.resources.snakefile.sequence_typing import get_sequence_typing_report, OUTPUT_TYPING_SUMMARY

#######################
# Included Snakefiles #
#######################
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: assembly_spades.SNAKEFILE_ASSEMBLY_SPADES
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING


rule all:
    """
    Rule to generate the required output files.
    """
    input:
        HTML = config['output_report'],
        TSV = config['output_tabular']

rule select_fastq:
    """
    This rule creates an IO object with the trimmed FASTQ files.
    Other workflows such as Kraken or Assembly rely on this dictionary to get input files (PE or SE).
    """
    input:
        FASTQ_PE = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_DICT,
    output:
        IO_FASTQ = Path(config['working_dir']) / 'fq_dict.io'
    shell:
        "cp {input.FASTQ_PE} {output.IO_FASTQ};"

rule select_fasta:
    """
    This rules links the output of the assembly workflow to the other workflows. 
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA
    output:
        FASTA = Path(config['working_dir']) / gene_detection.INPUT_GENE_DETECTION_FASTA
    shell:
        "cp {input.FASTA} {output.FASTA};"

rule report_pickle_citations:
    """
    This rule creates a pickle with a report section containing the citations.
    """
    output:
        IO = Path(config['working_dir']) / 'report' / 'html-citations.io'
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.components.html.htmlreportsection import HtmlReportSection
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        from camel.scripts.staphylococcuspipeline import CITATIONS_HTML
        section_citations = HtmlReportSection('Citations')
        with open(CITATIONS_HTML) as handle:
            section_citations.add_raw(handle.read())
        SnakemakeUtils.dump_object([ToolIOValue(section_citations)], output.IO)

rule report_command_section:
    """
    Creates a report section with the commands used in the pipeline. 
    """
    input:
        INFORMS_trimming = trimming.get_trimming_command_informs(config),
        INFORMS_assembly = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_INFORMS,
        INFORMS_assembly_filt = Path(config['working_dir']) / 'assembly_spades' / 'filtering' / 'informs.io',
        INFORMS_kraken = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS if 'kraken' in config['analyses'] else [],
        INFORMS_mapping = quality_checks.get_mapping_rate_informs(config),
        INFORMS_depth = quality_checks.get_depth_informs(config),
        # Gene detection
        INFORMS_argannot = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='argannot') if 'argannot' in config['analyses'] else [],
        INFORMS_card = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='card') if 'card' in config['analyses'] else [],
        INFORMS_ncbi_amr = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
        INFORMS_resfinder = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='resfinder') if 'resfinder' in config['analyses'] else [],
        INFORMS_virulencefinder = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        INFORMS_vfdb_core = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        INFORMS_plasmidfinder = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else [],
        # Sequence typing
        INFORMS_mlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        INFORMS_species = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='species_confirmation') if 'species_confirmation' in config['analyses'] else [],
        INFORMS_cgmlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        INFORMS_typing_amr = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='typing_amr') if 'typing_amr' in config['analyses'] else [],
        INFORMS_virulence = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='typing_virulence') if 'typing_virulence' in config['analyses'] else [],
        INFORMS_metal = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='metal_detergent') if 'metal_detergent' in config['analyses'] else [],
        INFORMS_pcr_sero = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='pcr_serogroup') if 'pcr_serogroup' in config['analyses'] else [],
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        working_dir = config['working_dir']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        informs = []
        for content in [SnakemakeUtils.load_object(io) for io in input]:
            if type(content) is dict:
                informs.append(content)
            elif type(content) is list:
                informs.extend(content)
        section = SnakePipelineUtils.create_commands_section(informs, params.working_dir)
        SnakemakeUtils.dump_object([ToolIOValue(section)], output.HTML)

rule report_combine_all:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        report_trimming = trimming.get_trimming_report(config),
        report_assembly = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_REPORT,
        report_kraken = Path(config['working_dir']) / (contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT if 'kraken' in config['analyses'] else contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY),
        report_adv_qc = Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT,
        # AMR detection
        report_argannot = gene_detection.get_gene_detection_report('argannot', config),
        report_card = gene_detection.get_gene_detection_report('card', config),
        report_ncbi_amr = gene_detection.get_gene_detection_report('ncbi_amr', config),
        report_resfinder = gene_detection.get_gene_detection_report('resfinder', config),
        # Virulence detection
        report_virulence = gene_detection.get_gene_detection_report('virulencefinder', config),
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        # Plasmid replicon detection
        report_plasmidfinder = gene_detection.get_gene_detection_report('plasmidfinder', config),
        # Typing
        report_mlst = get_sequence_typing_report('mlst', config),
        report_species = get_sequence_typing_report('species_confirmation', config),
        report_amr_typing = get_sequence_typing_report('typing_amr', config),
        report_cgmlst = get_sequence_typing_report('cgmlst', config),
        report_metal_detergent = get_sequence_typing_report('metal_detergent', config),
        report_pcr_serogroup = get_sequence_typing_report('pcr_serogroup', config),
        report_viru_typing = get_sequence_typing_report('typing_virulence', config),
        # Report
        report_citations = rules.report_pickle_citations.output.IO,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        fastq_input = config['fastq_pe'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        detection_method = config['detection_method']
    run:
        import datetime
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        # Add header section
        report = SnakePipelineUtils.init_pipeline_report(
            output.HTML, params.output_dir, params.pipeline_info)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            params.sample_name, datetime.datetime.now(), params.pipeline_info['version'],
            ', '.join(entry['name'] for entry in params.fastq_input), [('Detection method', params.detection_method)]))

        # Add report content
        report_structure = [
            ('Read trimming and basic QC', 'trim', [input.report_trimming]),
            ('Assembly', 'assem', [input.report_assembly]),
            ('Advanced QC', 'adv_qc', [input.report_kraken, input.report_adv_qc]),
            ('Species identification', 'species', [input.report_species, input.report_mlst]),
            ('AMR detection', 'amr', [
                input.report_argannot, input.report_card, input.report_ncbi_amr, input.report_resfinder]),
            ('Virulence detection', 'virulence', [input.report_virulence, input.report_vfdb_core]),
            ('Plasmid replicon detection', 'virulence', [input.report_plasmidfinder]),
            ('Sequence typing', 'typing', [
                input.report_amr_typing, input.report_cgmlst, input.report_metal_detergent,
                input.report_pcr_serogroup, input.report_viru_typing]),
            ('Citations', 'citations', [input.report_citations]),
            ('Commands', 'commands', [input.report_commands])
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
        input_filenames = ', '.join(entry['name'] for entry in config['fastq_pe' if 'fastq_pe' in config else 'fastq_se'])
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
    Combines the summary information of several steps into a single TSV file.
    """
    input:
        rules.summary_init.output.TSV,
        trimming.get_trimming_summary(config),
        Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_SUMMARY,
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_SUMMARY if 'kraken' in config['analyses'] else [],
        # AMR detection
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='argannot') if 'argannot' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='card') if 'card' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='resfinder') if 'resfinder' in config['analyses'] else [],
        # Virulence detection
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        # Plasmid replicon detection
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else [],
        # Sequence typing
        Path(config['working_dir']) / str(OUTPUT_TYPING_SUMMARY).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(OUTPUT_TYPING_SUMMARY).format(scheme='species_confirmation') if 'species_confirmation' in config['analyses'] else [],
        Path(config['working_dir']) / str(OUTPUT_TYPING_SUMMARY).format(scheme='typing_amr') if 'typing_amr' in config['analyses'] else [],
        Path(config['working_dir']) / str(OUTPUT_TYPING_SUMMARY).format(scheme='typing_virulence') if 'typing_amr' in config['analyses'] else [],
        Path(config['working_dir']) / str(OUTPUT_TYPING_SUMMARY).format(scheme='metal_detergent') if 'metal_detergent' in config['analyses'] else [],
        Path(config['working_dir']) / str(OUTPUT_TYPING_SUMMARY).format(scheme='pcr_serogroup') if 'pcr_serogroup' in config['analyses'] else [],
        Path(config['working_dir']) / str(OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
    output:
        TSV = config['output_tabular']
    run:
        with open(output[0], 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
