import re
from pathlib import Path

from camel.app.io.tooliovalue import ToolIOValue
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import trimming_illumina, trimming, quality_checks, deconseq
from camel.scripts.influenzapipeline.snakefile import genometyping_blastn, alignment, sequence_extraction

#######################
# Included Snakefiles #
#######################
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: deconseq.SNAKEFILE_DECONSEQ
include: genometyping_blastn.SNAKEFILE_GENOMETYPING
include: alignment.SNAKEFILE_MAPPING
include: sequence_extraction.SNAKEFILE_SEQ_EXTRACTION

##################
# AUXILIARY CODE #
##################
def get_genometyping_db_date(path_to_reference: str) -> str:
    """
    Retrieves the date the 'latest' symlink of the database points to
    :param path_to_reference:
    :return:
    """
    p = str(Path(path_to_reference).resolve())
    return re.search(r'(20[1-9][0-9][01][0-9][0-3][0-9])', p).groups()[0]

#########
# Rules #
#########
rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        config['output_report'],
         config['output_tabular'],
         Path(config['working_dir']) / 'fq_dict.io',
         Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_SAMTOOLS_DEPTH_ANALYZER_INFORMS,
         Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_ALIGNMENTSUMMARY,
         Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE_INDEX_PREFIX,
         Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE_ITERATIVE


rule select_fastq:
    """
    This rule creates an IO object with the trimmed FASTQ files.
    Other workflows such as Kraken or Assembly rely on this dictionary to get input files (PE or SE).
    """
    input:
        FASTQ_PE = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_DICT
    output:
        IO_FASTQ = Path(config['working_dir']) / 'fq_dict.io'
    shell:
        "cp {input.FASTQ_PE} {output.IO_FASTQ};"


rule report_command_section:
    input:
        INFORMS_trimming = trimming.get_trimming_command_informs(config),
         INFORMS_deconseq_pe_fwd = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_INFORMS_PE_FWD if 'deconseq' in config['analyses'] else [],
         INFORMS_deconseq_pe_rev = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_INFORMS_PE_REV if 'deconseq' in config['analyses'] else [],
         INFORMS_deconseq_se_fwd = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_INFORMS_SE_FWD if 'deconseq' in config['analyses'] else [],
         INFORMS_deconseq_se_rev = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_INFORMS_SE_REV if 'deconseq' in config['analyses'] else [],
         INFORMS_blast_genometyping =Path(config['working_dir']) / genometyping_blastn.OUTPUT_BLASTN_INFORMS if 'genometyping' in config['analyses'] else [],
         INFORMS_alignment =Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_INFORMS
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
        quality_checks = Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT,
        report_deconseq = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_REPORT if 'deconseq' in config['analyses'] else [],
        report_genometyping = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_REPORT if 'genometyping' in config['analyses'] else [],
        report_alignment = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_REPORT,
        # Report
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        fastq_input = config['fastq_pe'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline']
    run:
        import datetime

        # Add header section
        report = SnakePipelineUtils.init_pipeline_report(
            output.HTML, params.output_dir, params.pipeline_info)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            params.sample_name,
            datetime.datetime.now(),
            params.pipeline_info['version'], ', '.join(entry['name'] for entry in params.fastq_input),
            [('Database', config['genometyping_db_source'].upper()),
             ('Database date', get_genometyping_db_date(config['genometyping_db']))]))

        # Add output sections
        report_structure = [
            ('Read trimming and basic QC', 'trim', [input.report_trimming]),
            ('Advanced QC', 'qual', [input.quality_checks])
        ]
        if 'deconseq' in config['analyses']:
            report_structure.append(('Decontamination', 'deconseq', [input.report_deconseq]))
        if 'genometyping' in config['analyses']:
            report_structure.append(('Genome typing', 'genometyping', [input.report_genometyping]))
        report_structure.append(('Alignment', 'alignment', [input.report_alignment]))
        report_structure.append(('Commands', 'commands', [input.report_commands]))
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
                ('analysis_date', analysis_date)]:
                handle.write('\t'.join(kv_pair))
                handle.write('\n')

rule summary_combine_all:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        rules.summary_init.output.TSV,
        trimming.get_trimming_summary(config),
        Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_SUMMARY,
        Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_SUMMARY
    output:
        TSV = config['output_tabular']
    run:
        with open(output.TSV, 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
