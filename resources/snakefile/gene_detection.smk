import json
import shutil

from app.camel import Camel
from app.components.filesystemhelper import FileSystemHelper
from app.components.html.htmlreportsection import HtmlReportSection
from app.io.tooliodirectory import ToolIODirectory
from app.io.tooliovalue import ToolIOValue
from app.pipeline.snakestep import SnakeStep
from app.snakemake.snakemakeutils import SnakemakeUtils
from resources.snakefile.gene_detection import *

camel = Camel()


rule Gene_detection_db_manager:
    """
    Retrieves the FASTA file and the metadata from a database folder.
    """
    output:
        FASTA=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'fasta.io'),
        INFORMS=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'informs.io')
    params:
        db_path=lambda wildcards: config['gene_detection'][wildcards.db]['path'],
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}')
    run:
        from app.tools.pipelines.gene_detection.dbmanager import DBManager
        db_manager = DBManager(camel)
        db_manager.add_input_files({'DIR': [ToolIODirectory(params.db_path)]})
        SnakemakeUtils.add_pickle_inputs(db_manager, input)
        step = SnakeStep(rule, db_manager, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(db_manager, output)

rule Gene_detection_blastn:
    """
    Performs local alignment using Blastn+.
    """
    input:
        FASTA=os.path.join(config['working_dir'], INPUT_GENE_DETECTION_FASTA),
        DB_BLAST=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'fasta.io')
    output:
        ASN=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'blastn', 'asn.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'blastn')
    run:
        from app.tools.blast.blastn import Blastn
        blastn = Blastn(camel)
        SnakemakeUtils.add_pickle_inputs(blastn, input)
        step = SnakeStep(rule, blastn, camel, params.running_dir, config)
        blastn.update_parameters(threads=1)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blastn, output)

rule Gene_detection_tsv_generation:
    """
    Generates tabular output format to extract hit statistics.
    """
    input:
        ASN=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'blastn', 'asn.io')
    output:
        TSV=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'tsv_generation', 'tsv.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'tsv_generation')
    run:
        from app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter(camel)
        SnakemakeUtils.add_pickle_inputs(blast_formatter, input)
        step = SnakeStep(rule, blast_formatter, camel, params.running_dir, config)
        blast_formatter.update_parameters(output_format='"7 pident sseqid sseq slen qseqid qstart qend score"')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blast_formatter, output)

rule Gene_detection_hit_filtering:
    """
    Filters hits based on percent identity and query coverage.
    Extracts the hit information based on the database metadata.
    """
    input:
        TSV=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'tsv_generation', 'tsv.io'),
        INFORMS_db_info=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'informs.io')
    output:
        VAL_Hits=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_filtering', 'blast-hits.io'),
        TSV=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_filtering', 'tsv-filtered.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_filtering'),
        min_percent_identity=lambda wildcards: config['gene_detection'][wildcards.db].get('min_percent_identity'),
        min_coverage=lambda wildcards: config['gene_detection'][wildcards.db].get('min_coverage'),
        extra_column=lambda wildcards: config['gene_detection'][wildcards.db].get('extra_column'),
        output_filename = lambda wildcards: 'hits-{}-{}.tsv'.format(
            FileSystemHelper.make_valid(config['sample_name']),
            FileSystemHelper.make_valid(wildcards.db)),
        filtering_method=lambda wildcards: config['gene_detection'][wildcards.db].get('filtering_method'),
    run:
        from app.tools.pipelines.gene_detection.hitfiltering import HitFiltering
        hit_filtering = HitFiltering(camel)
        SnakemakeUtils.add_pickle_inputs(hit_filtering, input)
        step = SnakeStep(rule, hit_filtering, camel, params.running_dir, config)

        # Update parameters
        hit_filtering.update_parameters(output_filename=os.path.join(params.running_dir, params.output_filename))
        if params.filtering_method is not None:
            hit_filtering.update_parameters(filtering_method=params.filtering_method)
        if params.min_percent_identity is not None:
            hit_filtering.update_parameters(min_percent_identity=params.min_percent_identity)
        if params.min_coverage is not None:
            hit_filtering.update_parameters(min_coverage=params.min_coverage)
        if params.extra_column is not None:
           hit_filtering.update_parameters(
               extra_column_name=params.extra_column[0], extra_column_key=params.extra_column[1])

        # Run tool
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(hit_filtering, output)

rule Gene_detection_text_alignment_generation:
    """
    Generates alignments in the text format.
    """
    input:
        ASN=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'blastn', 'asn.io')
    output:
        TXT=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'alignment_generation', 'txt.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'alignment_generation')
    run:
        from app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter(camel)
        SnakemakeUtils.add_pickle_inputs(blast_formatter, input)
        step = SnakeStep(rule, blast_formatter, camel, params.running_dir, config)
        blast_formatter.update_parameters(output_format='0', num_alignments=1000)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blast_formatter, output)

rule Gene_detection_text_alignment_extraction:
    """
    Extracts a text alignment for the selected hits and attaches them to the hit objects.
    """
    input:
        TXT=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'alignment_generation', 'txt.io'),
        VAL_Hits=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_filtering', 'blast-hits.io')
    output:
        VAL_Hits=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'alignment_extraction', 'blast-hits.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'alignment_extraction')
    run:
        from app.tools.pipelines.gene_detection.alignmentextractor import AlignmentExtractor
        alignment_extractor = AlignmentExtractor(camel)
        SnakemakeUtils.add_pickle_inputs(alignment_extractor, input)
        step = SnakeStep(rule, alignment_extractor, camel, params.running_dir, config)
        step.run_step()
        hits_with_alignment = []
        for io_value, alignment in zip(SnakemakeUtils.load_object(input.VAL_Hits),
                                       alignment_extractor.tool_outputs['TXT']):
            io_value.value.alignment_path = alignment.path
            hits_with_alignment.append(io_value)
        SnakemakeUtils.dump_object(hits_with_alignment, output.VAL_Hits)

rule Gene_detection_srst2:
    """
    Read-mapping based gene detection using SRST2.
    """
    input:
        FASTQ_PE=os.path.join(config['working_dir'], INPUT_GENE_DETECTION_FASTQ_PE),
        FASTA=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'fasta.io')
    output:
        TSV=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'srst2', 'tsv.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'srst2')
    threads: 4
    run:
        from app.tools.srst2.srst2gene import Srst2Gene
        srst2 = Srst2Gene(camel)
        SnakemakeUtils.add_pickle_inputs(srst2, input)
        step = SnakeStep(rule, srst2, camel, params.running_dir, config)
        srst2.update_parameters(threads=threads, forward_designator='1P', reverse_designator='2P')
        step.run_step()
        if 'TSV' in srst2.tool_outputs:
            SnakemakeUtils.dump_tool_outputs(srst2, output)
        else:
            SnakemakeUtils.dump_object([], output.TSV)

rule Gene_detection_srst2_hit_extraction:
    """
    Extracts hits from the SRST2 output.
    """
    input:
        TSV=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'srst2', 'tsv.io'),
        INFORMS_db=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'informs.io')
    output:
        VAL_Hits=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'srst2', 'srst2-hits.io'),
        TSV=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_extraction', 'tsv-srst2.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'srst2'),
        extra_column=lambda wildcards: config['gene_detection'][wildcards.db].get('extra_column', None),
        output_filename = lambda wildcards: 'hits-{}-{}.tsv'.format(
            FileSystemHelper.make_valid(str(config['sample_name'])),
            FileSystemHelper.make_valid(wildcards.db))
    run:
        from app.tools.pipelines.gene_detection.srst2hitextractor import SRST2HitExtractor
        extractor = SRST2HitExtractor(camel)
        step = SnakeStep(rule, extractor, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(extractor, input)
        extractor.update_parameters(output_filename=os.path.join(params.running_dir, params.output_filename),
                                    extra_column=params.extra_column)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(extractor, output)

rule Gene_detection_get_hits:
    """
    Retrieves the hits from the blastn / SRST2 detection method based on the config
    """
    input:
        hits_blast=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'alignment_extraction', 'blast-hits.io') if config['detection_method'] == 'blast' else [],
        hits_srst2=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'srst2', 'srst2-hits.io') if config['detection_method'] == 'srst2' else [],
        tsv_blast=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_filtering', 'tsv-filtered.io') if config['detection_method'] == 'blast' else [],
        tsv_srst2=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_extraction', 'tsv-srst2.io') if config['detection_method'] == 'srst2' else []
    output:
        VAL_Hits=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_selection', 'selected-hits.io'),
        TSV=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_selection', 'selected-tsv.io')
    run:
        if len(input.hits_blast) > 0:
            shutil.copyfile(input.hits_blast, output.VAL_Hits)
        else:
            shutil.copyfile(input.hits_srst2, output.VAL_Hits)
        if len(input.tsv_blast) > 0:
            shutil.copyfile(input.tsv_blast, output.TSV)
        else:
            shutil.copyfile(input.tsv_srst2, output.TSV)

rule Gene_detection_report:
    """
    Creates HTML reports for the gene detection.
    """
    input:
        VAL_Hits=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_selection', 'selected-hits.io'),
        INFORMS_db_info=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'informs.io'),
        TSV=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_selection', 'selected-tsv.io')
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_REPORT),
    params:
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'report'),
        sample_name=config['sample_name']
    run:
        from app.tools.pipelines.gene_detection.htmlreportergenedetection import HtmlReporterGeneDetection
        reporter = HtmlReporterGeneDetection(camel)
        step = SnakeStep(rule, reporter, camel, params.running_dir, config)
        reporter.add_input_files({'SAMPLE_NAME': [ToolIOValue(params.sample_name)]})
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule Gene_detection_create_empty_report:
    """
    Creates an empty HTML report for the gene detection.
    """
    input:
        INFORMS_db_info=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'informs.io')
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_REPORT_EMPTY)
    params:
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'report'),
    run:
        db_info = SnakemakeUtils.load_object(input[0])
        section = HtmlReportSection(db_info['title'], 3)
        section.add_paragraph('Analysis disabled')
        SnakemakeUtils.dump_object([ToolIOValue(section)], output[0])

rule Gene_detection_dump_summary_info:
    """
    Dumps the summary information from the gene detection in tabular format.
    """
    input:
        INFORMS_hits=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'selected-hits.io')
    output:
        os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_SUMMARY)
    params:
        running_dir=os.path.join('gene_detection', '{db}', 'summary')
    run:
        informs = SnakemakeUtils.load_object(input.INFORMS_hits)
        hit_info = []
        for hit in informs:
            hit_info.append(hit.value.to_table_row().split('\t'))
        with open(output[0], 'w') as handle:
            handle.write('hits_{}\t{}'.format(wildcards.db, json.dumps(hit_info)))
            handle.write('\n')
