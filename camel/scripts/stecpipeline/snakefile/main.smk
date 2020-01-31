from pathlib import Path

from camel.resources.snakefile import trimming, trimming_illumina, trimming_iontorrent, assembly_spades, \
    quality_checks, contamination_check_kraken, gene_detection, pointfinder, variant_calling, variant_filtering, \
    sequence_typing
from camel.scripts.stecpipeline.snakefile import serotype_detection

#######################
# Included Snakefiles #
#######################
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: trimming_iontorrent.SNAKEFILE_TRIMMING_IONTORRENT
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: assembly_spades.SNAKEFILE_ASSEMBLY_SPADES
include: variant_calling.SNAKEFILE_VARIANT_CALLING
include: variant_filtering.SNAKEFILE_VARIANT_FILTERING
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: pointfinder.SNAKEFILE_POINTFINDER
include: serotype_detection.SNAKEFILE_SEROTYPE
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING


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


rule init_summary:
    """
    Initializes the summary output file.
    """
    output:
        summary = Path(config['working_dir']) / 'summary' / 'summary-init.tsv'
    run:
        import datetime
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        analysis_date = datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT)
        input_filenames = ', '.join(entry['name'] for entry in config['fastq_pe' if 'fastq_pe' in config else 'fastq_se'])
        with open(output.summary, 'w') as handle:
            for kv_pair in [
                ('pipeline_name', config['pipeline']['name']),
                ('pipeline_version', config['pipeline']['version']),
                ('sample', config['sample_name']),
                ('input_files', input_filenames),
                ('analysis_date', analysis_date),
                ('detection_method', config['detection_method']),
                ('read_type', config['read_type'])]:
                handle.write('\t'.join(kv_pair))
                handle.write('\n')


rule select_fastq:
    """
    This rule creates an IO object with the trimmed FASTQ files.
    Other workflows such as Kraken or Assembly rely on this dictionary to get input files (PE or SE).
    """
    input:
        FASTQ_PE = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_DICT if trimming.get_read_type(config) == 'illumina' else [],
        FASTQ_SE = Path(config['working_dir']) / trimming_iontorrent.OUTPUT_TRIMMING_IONTORRENT_DICT if trimming.get_read_type(config) == 'iontorrent' else []
    output:
        IO_FASTQ = Path(config['working_dir']) / 'fq_dict.io'
    params:
        read_type = config['read_type']
    run:
        import shutil
        for key, fq in input.items():
            if len(fq) == 0:
                continue
            shutil.copyfile(fq, output.IO_FASTQ)


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
        from camel.scripts.stecpipeline import CITATIONS_HTML
        section_citations = HtmlReportSection('Citations')
        with open(CITATIONS_HTML) as handle:
            section_citations.add_raw(handle.read())
        SnakemakeUtils.dump_object([ToolIOValue(section_citations)], output.IO)


rule report_command_section:
    input:
        INFORMS_trimming = trimming.get_trimming_command_informs(config),
        INFORMS_assembly = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_INFORMS,
        INFORMS_assembly_filt=Path(config['working_dir']) / 'assembly_spades' / 'filtering' / 'informs.io',
        INFORMS_kraken = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS if 'kraken' in config['analyses'] else [],
        INFORMS_mapping = quality_checks.get_mapping_rate_informs(config),
        INFORMS_depth = quality_checks.get_depth_informs(config),
        INFORMS_variant_calling_all = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_INFORMS_ALL,
        INFORMS_variant_filtering_all = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_INFORMS_ALL,
        INFORMS_resfinder = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='resfinder') if 'resfinder' in config['analyses'] else [],
        INFORMS_argannot = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='argannot') if 'argannot' in config['analyses'] else [],
        INFORMS_card = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='card') if 'card' in config['analyses'] else [],
        INFORMS_ncbi_amr = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
        INFORMS_pointfinder = Path(config['working_dir']) / pointfinder.OUTPUT_POINTFINDER_INFORMS if 'pointfinder' in config['analyses'] else [],
        INFORMS_virulence = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        INFORMS_virulence_shiga = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='virulencefinder_shiga') if 'virulencefinder' in config['analyses'] else [],
        INFORMS_serotype_h = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='serotype_h') if 'serotype' in config['analyses'] else [],
        INFORMS_serotype_o = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='serotype_o') if 'serotype' in config['analyses'] else [],
        INFORMS_plasmidfinder = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else []
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        working_dir = config['working_dir']
    run:
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


rule combine_reports:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        report_trimming=trimming.get_trimming_report(config),
        report_assembly=Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_REPORT,
        report_kraken=Path(config['working_dir']) / (contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT if 'kraken' in config['analyses'] else contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY),
        report_adv_qc=Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT,
        report_variant=Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_REPORT,
        report_pointfinder=Path(config['working_dir']) / (pointfinder.OUTPUT_POINTFINDER_REPORT if 'pointfinder' in config['analyses'] else pointfinder.OUTPUT_POINTFINDER_REPORT_EMPTY),
        report_serotype=Path(config['working_dir']) / serotype_detection.OUTPUT_SEROTYPE_REPORT,
        # Gene detection
        report_resfinder=gene_detection.get_gene_detection_report('resfinder', config),
        report_argannot=gene_detection.get_gene_detection_report('argannot', config),
        report_card=gene_detection.get_gene_detection_report('card', config),
        report_ncbi_amr=gene_detection.get_gene_detection_report('ncbi_amr', config),
        report_virulence=gene_detection.get_gene_detection_report('virulencefinder', config),
        report_virulence_shiga=gene_detection.get_gene_detection_report('virulencefinder_shiga', config, 'virulencefinder'),
        report_plasmidfinder=gene_detection.get_gene_detection_report('plasmidfinder', config),
        report_serotype_o_type=gene_detection.get_gene_detection_report('serotype_o', config, 'serotype'),
        report_serotype_h_type=gene_detection.get_gene_detection_report('serotype_h', config, 'serotype'),
        # Typing
        report_mlst_warwick=sequence_typing.get_sequence_typing_report('mlst_warwick', config),
        report_mlst_pasteur=sequence_typing.get_sequence_typing_report('mlst_pasteur', config),
        report_cgmlst=sequence_typing.get_sequence_typing_report('cgmlst', config),
        # Report
        report_citations = rules.report_pickle_citations.output.IO,
        report_commands = rules.report_command_section.output.HTML
    output:
        report = config['output_report']
    params:
        sample_name = config['sample_name'],
        fastq_input = config['fastq_pe' if 'fastq_pe' in config else 'fastq_se'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        detection_method = config['detection_method'],
        read_type = config['read_type']
    run:
        import datetime
        from camel.app.components.html.htmlreport import HtmlReport
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.resources import CSS_STYLE
        from camel.resources.javascript import JQUERY_SRC

        # Add header section
        report = HtmlReport(output.report, params.output_dir, [JQUERY_SRC])
        report.initialize(params.pipeline_info['name'], CSS_STYLE)
        report.add_pipeline_header(f"{params.pipeline_info['name']} {params.pipeline_info['version']}")
        report.add_html_object(SnakePipelineUtils.create_input_section(
            params.sample_name,
            datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT),
            params.pipeline_info['version'],
            ', '.join(entry['name'] for entry in params.fastq_input),
            [('Detection method', params.detection_method), ('Read type', params.read_type)],
        ))

        # Add output sections
        report_structure = [
            ('Read trimming and basic QC', 'trim', [input.report_trimming]),
            ('Assembly', 'assem', [input.report_assembly]),
            ('Advanced QC', 'adv_qc', [input.report_kraken, input.report_adv_qc]),
            ('Variant calling', 'variant', [input.report_variant]),
            ('Resistance characterization', 'res', [input.report_resfinder, input.report_argannot, input.report_card,
                                                    input.report_ncbi_amr, input.report_pointfinder]),
            ('Virulence characterization', 'viru', [input.report_virulence, input.report_virulence_shiga]),
            ('Serotype determination', 'sero', [input.report_serotype_o_type, input.report_serotype_h_type,
                                                input.report_serotype]),
            ('Plasmid replicon detection', 'plasmid', [input.report_plasmidfinder]),
            ('Sequence typing', 'st', [input.report_mlst_warwick, input.report_mlst_pasteur, input.report_cgmlst]),
            ('Citations', 'citations', [input.report_citations]),
            ('Commands', 'commands', [input.report_commands])
        ]
        SnakePipelineUtils.add_report_content(report, report_structure)


rule combine_summary_files:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        rules.init_summary.output.TSV,
        trimming.get_trimming_summary(config),
        Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_SUMMARY,
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_SUMMARY if 'kraken' in config['analyses'] else [],
        Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_SUMMARY,
        Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_SUMMARY,
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='resfinder') if 'resfinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='card') if 'card' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='argannot') if 'argannot' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
        Path(config['working_dir']) / pointfinder.OUTPUT_POINTFINDER_SUMMARY if 'pointfinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='virulencefinder_shiga') if 'virulencefinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='serotype_h') if 'serotype' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='serotype_o') if 'serotype' in config['analyses'] else [],
        Path(config['working_dir']) / serotype_detection.OUTPUT_SEROTYPE_SUMMARY if 'serotype' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst_pasteur') if 'mlst_pasteur' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst_warwick') if 'mlst_warwick' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else []
    output:
        config.get('output_tabular')
    run:
        with open(output[0], 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
