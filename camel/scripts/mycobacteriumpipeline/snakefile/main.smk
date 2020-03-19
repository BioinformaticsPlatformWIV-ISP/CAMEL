from pathlib import Path

from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import trimming_illumina, contamination_check_kraken, quality_checks, assembly_spades, \
    variant_calling, variant_filtering, gene_detection, sequence_typing, trimming, pointfinder
from camel.scripts.mycobacteriumpipeline.snakefile import csb_rd, snpit, hsp65, spoligotyping, snplineage, assay51snp, \
    amr

#######################
# Included Snakefiles #
#######################
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: assembly_spades.SNAKEFILE_ASSEMBLY_SPADES
include: variant_calling.SNAKEFILE_VARIANT_CALLING
include: variant_filtering.SNAKEFILE_VARIANT_FILTERING
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING
include: csb_rd.SNAKEFILE_CSB_RD
include: snpit.SNAKEFILE_SNPIT
include: hsp65.SNAKEFILE_HSP65
include: assay51snp.SNAKEFILE_51SNP
include: spoligotyping.SNAKEFILE_SPOLIGOTYPING
include: snplineage.SNAKEFILE_SNP_LINEAGE
include: amr.SNAKEFILE_AMR
include: pointfinder.SNAKEFILE_POINTFINDER


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
        from camel.scripts.mycobacteriumpipeline import CITATIONS_HTML
        section_citations = HtmlReportSection('Citations')
        with open(CITATIONS_HTML) as handle:
            section_citations.add_raw(handle.read())
        SnakemakeUtils.dump_object([ToolIOValue(section_citations)], output[0])

rule report_command_section:
    input:
        INFORMS_trimming = trimming.get_trimming_command_informs(config),
        INFORMS_assembly = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_INFORMS,
        INFORMS_assembly_filt = Path(config['working_dir']) / 'assembly_spades' / 'filtering' / 'informs.io',
        INFORMS_kraken = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS if 'kraken' in config['analyses'] else [],
        INFORMS_mapping = quality_checks.get_mapping_rate_informs(config),
        INFORMS_depth = quality_checks.get_depth_informs(config),
        INFORMS_variant_calling_all = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_INFORMS_ALL,
        INFORMS_variant_filtering_all = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_INFORMS_ALL,
        INFORMS_snpit = Path(config['working_dir']) / snpit.OUTPUT_SNPIT_INFORMS if 'snpit' in config['analyses'] else [],
        INFORMS_16s = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='ncbi_16s') if 'ncbi_16s' in config['analyses'] else [],
        INFORMS_csb_rd = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='csb_rd') if 'csb_rd,' in config['analyses'] else [],
        INFORMS_hsp65 = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='hsp65') if 'hsp65' in config['analyses'] else [],
        INFORMS_spoligo = Path(config['working_dir']) / spoligotyping.OUTPUT_SPOLIGOTYPING_INFORMS if 'spoligotyping' in config['analyses'] else [],
        INFORMS_pointfinder = Path(config['working_dir']) / pointfinder.OUTPUT_POINTFINDER_INFORMS if 'pointfinder' in config['analyses'] else []
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        working_dir = config['working_dir']
    run:
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
        report_variant = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_REPORT,
        report_ncbi_16s = gene_detection.get_gene_detection_report('ncbi_16s', config),
        report_51snp = Path(config['working_dir']) / (assay51snp.OUTPUT_51SNP_REPORT if '51snp' in config['analyses'] else assay51snp.OUTPUT_51SNP_REPORT_EMPTY),
        report_csb_rd = Path(config['working_dir']) / (csb_rd.OUTPUT_CSB_RD_REPORT if 'csb_rd' in config['analyses'] else csb_rd.OUTPUT_CSB_RD_REPORT_EMPTY),
        report_hsp65 = Path(config['working_dir']) / (hsp65.OUTPUT_HSP65_REPORT if 'hsp65' in config['analyses'] else hsp65.OUTPUT_HSP65_REPORT_EMPTY),
        report_snpit = Path(config['working_dir']) / (snpit.OUTPUT_SNPIT_REPORT if 'snpit' in config['analyses'] else snpit.OUTPUT_SNPIT_REPORT_EMPTY),
        report_spoligo = Path(config['working_dir']) / (spoligotyping.OUTPUT_SPOLIGOTYPING_REPORT if 'spoligotyping' in config['analyses'] else spoligotyping.OUTPUT_SPOLIGOTYPING_REPORT_EMPTY),
        report_snp_lineage = Path(config['working_dir']) / (snplineage.OUTPUT_SNP_LINEAGE_REPORT if 'snp_lineage' in config['analyses'] else snplineage.OUTPUT_SNP_LINEAGE_REPORT_EMPTY),
        report_amr = Path(config['working_dir']) / (amr.OUTPUT_AMR_REPORT if 'amr' in config['analyses'] else amr.OUTPUT_AMR_REPORT_EMPTY),
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        # Report
        report_citations = rules.report_pickle_citations.output.IO,
        report_commands = rules.report_command_section.output.HTML,
        report_pointfinder = Path(config['working_dir']) / (pointfinder.OUTPUT_POINTFINDER_REPORT if 'pointfinder' in config['analyses'] else pointfinder.OUTPUT_POINTFINDER_REPORT_EMPTY),
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

        # Add header section
        report = SnakePipelineUtils.init_pipeline_report(
            output.HTML, params.output_dir, params.pipeline_info)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            params.sample_name,
            datetime.datetime.now(),
            params.pipeline_info['version'], ', '.join(entry['name'] for entry in params.fastq_input),
            [('Detection method', params.detection_method)]))

        # Add output sections
        report_structure = [
            ('Read trimming and basic QC', 'trim', [input.report_trimming]),
            ('Assembly', 'assembly', [input.report_assembly]),
            ('Advanced QC', 'adv_qc', [input.report_kraken, input.report_adv_qc]),
            ('Variant calling', 'variant', [input.report_variant]),
            ('Species identification', 'identification', [
                input.report_ncbi_16s, input.report_snpit, input.report_csb_rd, input.report_hsp65,
                input.report_51snp]),
            ('PointFinder', 'pointfinder', [input.report_pointfinder]),
            ('Spoligotyping and lineage', 'spoligotyping', [input.report_spoligo, input.report_snp_lineage]),
            ('AMR detection', 'amr', [input.report_amr]),
            ('Sequence typing', 'typing', [input.report_mlst, input.report_cgmlst]),
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
        analysis_date = datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT)
        input_filenames = ', '.join(entry['name'] for entry in config['fastq_pe'])
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
        trimming.get_trimming_summary(config),
        Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_SUMMARY,
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_SUMMARY if 'kraken' in config['analyses'] else [],
        Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_SUMMARY,
        Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_SUMMARY,
        Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_SUMMARY if 'kraken' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='ncbi_16s') if 'ncbi_16s' in config['analyses'] else [],
        Path(config['working_dir']) / csb_rd.OUTPUT_CSB_RD_SUMMARY if 'csb_rd' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='hsp65') if 'hsp65' in config['analyses'] else [],
        Path(config['working_dir']) / assay51snp.OUTPUT_51SNP_SUMMARY if '51snp' in config['analyses'] else [],
        Path(config['working_dir']) / snpit.OUTPUT_SNPIT_SUMMARY if 'snpit' in config['analyses'] else [],
        Path(config['working_dir']) / spoligotyping.OUTPUT_SPOLIGOTYPING_SUMMARY if 'spoligotyping' in config['analyses'] else [],
        Path(config['working_dir']) / snplineage.OUTPUT_SNP_LINEAGE_SUMMARY if 'snp_lineage' in config['analyses'] else [],
        Path(config['working_dir']) / amr.OUTPUT_AMR_SUMMARY if 'amr' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        Path(config['working_dir']) / pointfinder.OUTPUT_POINTFINDER_SUMMARY if 'pointfinder' in config['analyses'] else []
    output:
        TSV = config['output_tabular']
    run:
        with open(output.TSV, 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
