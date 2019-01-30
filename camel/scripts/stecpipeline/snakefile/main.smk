import os

from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

from camel.scripts.stecpipeline import SNAKEFILE_SEROTYPE
from camel.resources.snakefile import SNAKEFILE_READ_TRIMMING, SNAKEFILE_ASSEMBLY_SPADES, SNAKEFILE_GENE_DETECTION, \
    SNAKEFILE_SEQUENCE_TYPING, SNAKEFILE_VARIANT_CALLING, SNAKEFILE_VARIANT_FILTERING, \
    SNAKEFILE_CONTAMINATION_CHECK_KRAKEN, SNAKEFILE_POINTFINDER, SNAKEFILE_ADV_QC
from camel.resources.snakefile.assembly_spades import OUTPUT_ASSEMBLY_REPORT, OUTPUT_ASSEMBLY_FASTA
from camel.resources.snakefile.contamination_check_kraken import OUTPUT_CONTAMINATION_CHECK_REPORT, \
    OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY
from camel.resources.snakefile.gene_detection import INPUT_GENE_DETECTION_FASTA, get_gene_detection_report
from camel.resources.snakefile.pointfinder import OUTPUT_POINTFINDER_REPORT, OUTPUT_POINTFINDER_REPORT_EMPTY
from camel.resources.snakefile.quality_checks import OUTPUT_QUALITY_CHECKS_REPORT
from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_REPORT
from camel.resources.snakefile.sequence_typing import get_sequence_typing_report
from camel.resources.snakefile.variant_calling import OUTPUT_VARIANT_CALLING_REPORT

#######################
# Included Snakefiles #
#######################
from camel.scripts.stecpipeline.snakefile.serotype_detection import OUTPUT_SEROTYPE_REPORT

include: SNAKEFILE_READ_TRIMMING
include: SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: SNAKEFILE_ADV_QC
include: SNAKEFILE_ASSEMBLY_SPADES
include: SNAKEFILE_VARIANT_CALLING
include: SNAKEFILE_VARIANT_FILTERING
include: SNAKEFILE_GENE_DETECTION
include: SNAKEFILE_POINTFINDER
include: SNAKEFILE_SEROTYPE
include: SNAKEFILE_SEQUENCE_TYPING


#########
# Rules #
#########
rule All:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        config['output_report'],
        config['output_tabular']


rule Init_summary:
    """
    Initializes the summary output file.
    """
    output:
        summary = os.path.join(config['working_dir'], 'summary', 'summary-init.tsv')
    run:
        import datetime
        analysis_date = datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT)
        input_filenames = ', '.join(entry['name'] for entry in config['fastq_pe'])
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


rule Combine_reports:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        report_trimming=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_REPORT),
        report_assembly=os.path.join(config['working_dir'], OUTPUT_ASSEMBLY_REPORT),
        report_kraken=os.path.join(config['working_dir'], OUTPUT_CONTAMINATION_CHECK_REPORT) if 'kraken' in config['analyses'] else os.path.join(config['working_dir'], OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY),
        report_adv_qc=os.path.join(config['working_dir'], OUTPUT_QUALITY_CHECKS_REPORT),
        report_variant=os.path.join(config['working_dir'], OUTPUT_VARIANT_CALLING_REPORT),
        report_resfinder=get_gene_detection_report('resfinder', config),
        report_argannot=get_gene_detection_report('argannot', config),
        report_card=get_gene_detection_report('card', config),
        report_ncbi_amr=get_gene_detection_report('ncbi_amr', config),
        report_pointfinder=os.path.join(config['working_dir'], OUTPUT_POINTFINDER_REPORT) if 'pointfinder' in config['analyses'] else os.path.join(config['working_dir'], OUTPUT_POINTFINDER_REPORT_EMPTY),
        report_virulence=get_gene_detection_report('virulencefinder', config),
        report_virulence_shiga=get_gene_detection_report('virulencefinder_shiga', config),
        report_plasmidfinder=get_gene_detection_report('plasmidfinder', config),
        report_serotype_o_type=get_gene_detection_report('serotype_o', config, 'serotype'),
        report_serotype_h_type=get_gene_detection_report('serotype_h', config, 'serotype'),
        report_serotype=os.path.join(config['working_dir'], OUTPUT_SEROTYPE_REPORT),
        report_mlst_warwick=get_sequence_typing_report('mlst_warwick', config),
        report_mlst_pasteur=get_sequence_typing_report('mlst_pasteur', config),
        report_cgmlst=get_sequence_typing_report('cgmlst', config)
    output:
        report = config['output_report']
    params:
        sample_name = config['sample_name'],
        fastq_input = config['fastq_pe'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        detection_method = config['detection_method'],
        read_type = config['read_type']
    run:
        import datetime
        from camel.app.components.html.htmlelement import HtmlElement
        from camel.app.components.html.htmlreport import HtmlReport
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
            ('Sequence typing', 'st', [input.report_mlst_warwick, input.report_mlst_pasteur, input.report_cgmlst])
        ]

        report.add_module_header('Sections')
        overview_list = HtmlElement('ul')
        for title, key, _ in report_structure:
            list_item = HtmlElement('li')
            list_item.add_html_object(HtmlElement('a', title, [('href', '#{}'.format(key))]))
            overview_list.add_html_object(list_item)
        report.add_html_object(overview_list)

        for title, key, items in report_structure:
            report.add_module_header(title, key)
            for pickle in items:
                if len(pickle) == 0:
                    continue
                section = SnakemakeUtils.load_object(pickle)[0].value
                report.add_html_object(section)
                section.copy_files(params.output_dir)

        report.save()


rule Combine_summary_files:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        os.path.join(config['working_dir'], 'summary', 'summary-init.tsv'),
    output:
        config.get('output_tabular')
    run:
        with open(output[0], 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())


rule Link_assembly_gene_detection:
    """
    Links the output of the assembly to the input of the gene detection.
    """
    input:
        os.path.join(config['working_dir'], OUTPUT_ASSEMBLY_FASTA)
    output:
        os.path.join(config['working_dir'], INPUT_GENE_DETECTION_FASTA)
    shell:
        "cp {input} {output}"
