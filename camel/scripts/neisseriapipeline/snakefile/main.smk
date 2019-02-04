import os

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import SNAKEFILE_READ_TRIMMING, SNAKEFILE_ASSEMBLY_SPADES, SNAKEFILE_GENE_DETECTION, \
    SNAKEFILE_SEQUENCE_TYPING, SNAKEFILE_CONTAMINATION_CHECK_KRAKEN, SNAKEFILE_ADV_QC
from camel.resources.snakefile.assembly_spades import OUTPUT_ASSEMBLY_REPORT, OUTPUT_ASSEMBLY_FASTA, \
    OUTPUT_ASSEMBLY_SUMMARY
from camel.resources.snakefile.contamination_check_kraken import OUTPUT_CONTAMINATION_CHECK_REPORT, \
    OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY, OUTPUT_CONTAMINATION_SUMMARY
from camel.scripts.neisseriapipeline import SNAKEFILE_SEROGROUP_DETERMINATION
from camel.scripts.neisseriapipeline.snakefile.serogroup_determination import OUTPUT_SEROGROUP_DETERMINATION_REPORT, \
    OUTPUT_SEROGROUP_DETERMINATION_REPORT_EMPTY, OUTPUT_SEROGROUP_DETERMINATION_SUMMARY
from camel.resources.snakefile.gene_detection import INPUT_GENE_DETECTION_FASTA, get_gene_detection_report, \
    OUTPUT_GENE_DETECTION_SUMMARY, INPUT_GENE_DETECTION_FASTQ
from camel.resources.snakefile.quality_checks import OUTPUT_QUALITY_CHECKS_REPORT, OUTPUT_QUALITY_CHECKS_SUMMARY
from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_REPORT, OUTPUT_READ_TRIMMING_SUMMARY, \
    OUTPUT_READ_TRIMMING_READS_PE
from camel.resources.snakefile.sequence_typing import get_sequence_typing_report, OUTPUT_TYPING_SUMMARY

#######################
# Included Snakefiles #
#######################
include: SNAKEFILE_READ_TRIMMING
include: SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: SNAKEFILE_ADV_QC
include: SNAKEFILE_ASSEMBLY_SPADES
include: SNAKEFILE_GENE_DETECTION
include: SNAKEFILE_SEQUENCE_TYPING
include: SNAKEFILE_SEROGROUP_DETERMINATION


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
                ('detection_method', config['detection_method'])]:
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
        report_resfinder=get_gene_detection_report('resfinder', config),
        report_argannot=get_gene_detection_report('argannot', config),
        report_card=get_gene_detection_report('card', config),
        report_ncbi_amr=get_gene_detection_report('ncbi_amr', config),
        report_mlst=get_sequence_typing_report('mlst', config),
        report_cgmlst=get_sequence_typing_report('cgmlst', config),
        report_pora=get_sequence_typing_report('pora', config),
        report_porb=get_sequence_typing_report('porb', config),
        report_feta=get_sequence_typing_report('feta', config),
        report_rplf=get_sequence_typing_report('rplf', config),
        report_vaccine_targets=get_sequence_typing_report('vaccine_targets', config),
        report_resistance_genes=os.path.join(config['working_dir'], 'typing', 'resistance_genes', 'metadata', 'report.html') if 'resistance_genes' in config['analyses'] else os.path.join(config['working_dir'], OUTPUT_TYPING_REPORT_EMPTY.format(scheme='resistance_genes')),
        report_fhbp=get_sequence_typing_report('fhbp', config),
        report_bast=get_sequence_typing_report('bast', config),
        report_serogroup=os.path.join(config['working_dir'], OUTPUT_SEROGROUP_DETERMINATION_REPORT if 'serogroup' in config['analyses'] else OUTPUT_SEROGROUP_DETERMINATION_REPORT_EMPTY),
        report_citations=os.path.join(config['working_dir'], 'report', 'html-citations.io')
    output:
        report = config['output_report']
    params:
        sample_name = config['sample_name'],
        fastq_input = config['fastq_pe'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        detection_method = config['detection_method'],
    run:
        import datetime
        from camel.app.components.html.htmlelement import HtmlElement
        from camel.app.components.html.htmlreport import HtmlReport
        from camel.resources import CSS_STYLE
        from camel.resources.javascript import JQUERY_SRC

        # Add header section
        report = HtmlReport(output.report, params.output_dir, [JQUERY_SRC])
        report.initialize(params.pipeline_info['name'], CSS_STYLE)
        report.add_pipeline_header(f"{params.pipeline_info['title']} {params.pipeline_info['version']}")
        report.add_html_object(SnakePipelineUtils.create_input_section(
            params.sample_name,
            datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT),
            params.pipeline_info['version'], ', '.join(entry['name'] for entry in params.fastq_input),
            [('Detection method', params.detection_method)]))

        # Add output sections
        report_structure = [
            ('Read trimming and basic QC', 'trim', [input.report_trimming]),
            ('Assembly', 'assem', [input.report_assembly]),
            ('Advanced QC', 'adv_qc', [input.report_kraken, input.report_adv_qc]),
            ('Resistance characterization', 'res', [input.report_resfinder, input.report_argannot, input.report_card,
                                                    input.report_ncbi_amr]),
            ('Sequence typing', 'st', [input.report_mlst, input.report_rplf, input.report_bast, input.report_pora,
                                       input.report_porb, input.report_feta, input.report_resistance_genes,
                                       input.report_vaccine_targets, input.report_fhbp, input.report_cgmlst]),
            ('Serogroup determination', 'serogroup', [input.report_serogroup]),
            ('Citations', 'citations', [input.report_citations])
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


rule Pickle_citations:
    """
    This rule creates a pickle with a report section containing the citations.
    """
    output:
        os.path.join(config['working_dir'], 'report', 'html-citations.io')
    run:
        from camel.scripts.neisseriapipeline import CITATIONS_HTML
        section_citations = HtmlReportSection('Citations')
        with open(CITATIONS_HTML) as handle:
            section_citations.add_raw(handle.read())
        SnakemakeUtils.dump_object([ToolIOValue(section_citations)], output[0])


rule Combine_summary_files:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        os.path.join(config['working_dir'], 'summary', 'summary-init.tsv'),
        os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_SUMMARY),
        os.path.join(config['working_dir'], OUTPUT_ASSEMBLY_SUMMARY),
        os.path.join(config['working_dir'], OUTPUT_QUALITY_CHECKS_SUMMARY),
        os.path.join(config['working_dir'], OUTPUT_CONTAMINATION_SUMMARY) if 'kraken' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_SUMMARY.format(db='resfinder')) if 'resfinder' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_SUMMARY.format(db='card')) if 'card' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_SUMMARY.format(db='argannot')) if 'argannot' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_SUMMARY.format(db='ncbi_amr')) if 'ncbi_amr' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_TYPING_SUMMARY.format(scheme='mlst')) if 'mlst' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_TYPING_SUMMARY.format(scheme='rplf')) if 'rplf' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_TYPING_SUMMARY.format(scheme='bast')) if 'bast' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_TYPING_SUMMARY.format(scheme='pora')) if 'pora' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_TYPING_SUMMARY.format(scheme='porb')) if 'porb' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_TYPING_SUMMARY.format(scheme='feta')) if 'feta' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_TYPING_SUMMARY.format(scheme='fhbp')) if 'fhbp' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_TYPING_SUMMARY.format(scheme='resistance_genes')) if 'resistance_genes' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_TYPING_SUMMARY.format(scheme='vaccine_targets')) if 'vaccine_targets' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_TYPING_SUMMARY.format(scheme='cgmlst')) if 'cgmlst' in config['analyses'] else [],
        os.path.join(config['working_dir'], OUTPUT_SEROGROUP_DETERMINATION_SUMMARY) if 'serogroup' in config['analyses'] else []
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

rule Link_trimming_gene_detection:
    """
    Links the output of the assembly to the input of the gene detection.
    """
    input:
        os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_PE)
    output:
        os.path.join(config['working_dir'], INPUT_GENE_DETECTION_FASTQ)
    shell:
        "cp {input} {output}"

rule Collect_read_mapping_input:
    """
    Collects the input for the read mapping.
    """
    input:
        ILLUMINA_FASTQ_PE = os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_PE) if config.get('read_type', 'illumina') == 'illumina' else [],
        ILLUMINA_FASTQ_SE_FWD = os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_SE_FWD) if config.get('read_type', 'illumina') == 'illumina' else [],
        ILLUMINA_FASTQ_SE_REV = os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_SE_REV) if config.get('read_type', 'illumina') == 'illumina' else []
    output:
        FASTQ=os.path.join(config['working_dir'], 'variant_calling', 'input-fastq.io')
    params:
        read_type = config.get('read_type', 'illumina')
    run:
        output_dict = {}
        if params.read_type == 'illumina':
            output_dict = {'FASTQ_PE': SnakemakeUtils.load_object(input.ILLUMINA_FASTQ_PE)}
            se_reads = SnakemakeUtils.load_object(input.ILLUMINA_FASTQ_SE_FWD) + \
                       SnakemakeUtils.load_object(input.ILLUMINA_FASTQ_SE_REV)
            if len(se_reads) > 0:
                output_dict['FASTQ_SE'] = se_reads
        else:
            output_dict = {'FASTQ_SE': SnakemakeUtils.load_object(input.IONTORRENT_FASTQ_SE)}
        SnakemakeUtils.dump_object(output_dict, output[0])


rule Parse_additional_resistance_gene_metadata:
    """
    This rule is used to add resistance gene metadata for penA and rpoB genes.
    The data is parsed from the PubMLST webpage.
    """
    input:
        hits=os.path.join(config['working_dir'], OUTPUT_TYPING_HITS.format(scheme='resistance_genes', locus_type='DNA', detection_method=config['detection_method'])),
        VAL_HTML=get_sequence_typing_report('resistance_genes', config),
        INFORMS_scheme = os.path.join(config['working_dir'], 'typing', 'resistance_genes', 'informs-locus_set.io')
    output:
        VAL_HTML=os.path.join(config['working_dir'], 'typing', 'resistance_genes', 'metadata', 'report.html')
    params:
        working_dir=os.path.join(config['working_dir'], 'typing', 'resistance_genes', 'metadata'),
        loci='penA, rpoB'
    run:
        from camel.app.tools.pipelines.neisseria.resistancemetadataextractor import ResistanceMetadataExtractor
        extractor = ResistanceMetadataExtractor(camel)
        SnakemakeUtils.add_pickle_inputs(extractor, input)
        step = Step(rule, extractor, camel, params.working_dir, config)
        extractor.update_parameters(loci=params.loci)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(extractor, output)
