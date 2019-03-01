import json
import shutil

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.genedetection.genedetectionutils import GeneDetectionUtils
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.gene_detection import *

camel = Camel.get_instance()


rule Gene_detection_db_manager:
    """
    Retrieves the FASTA file and the metadata from a database folder.
    """
    output:
        FASTA=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'fasta.io'),
        FASTA_clustered=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'fasta-clust.io'),
        INFORMS=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'informs.io')
    params:
        db_path=lambda wildcards: config['gene_detection'][wildcards.db]['path'],
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}')
    run:
        from camel.app.tools.pipelines.genedetection.dbmanager import DBManager
        db_manager = DBManager(camel)
        db_manager.add_input_files({'DIR': [ToolIODirectory(params.db_path)]})
        step = Step(rule, db_manager, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(db_manager, output)

rule Gene_detection_blastn:
    """
    Performs local alignment using Blastn+.
    """
    input:
        FASTA=os.path.join(config['working_dir'], INPUT_GENE_DETECTION_FASTA),
        DB_BLAST=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'fasta-clust.io')
    output:
        ASN = os.path.join(config['working_dir'], 'gene_detection', '{db}', 'blastn', 'asn.io'),
        INFORMS = os.path.join(config['working_dir'], 'gene_detection', '{db}', 'blastn', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'blastn')
    run:
        from camel.app.tools.blast.blastn import Blastn
        blastn = Blastn(camel)
        SnakemakeUtils.add_pickle_inputs(blastn, input)
        step = Step(rule, blastn, camel, params.running_dir, config)
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
        from camel.app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter(camel)
        SnakemakeUtils.add_pickle_inputs(blast_formatter, input)
        step = Step(rule, blast_formatter, camel, params.running_dir, config)
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
        TSV=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_filtering', 'tsv-filtered.io'),
        INFORMS=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_filtering', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_filtering'),
        output_filename = lambda wildcards: 'hits-{}-{}.tsv'.format(
            FileSystemHelper.make_valid(config['sample_name']), FileSystemHelper.make_valid(wildcards.db)),
        db_config = lambda wildcards: config['gene_detection'][wildcards.db]
    run:
        from camel.app.tools.pipelines.genedetection.blasthitfiltering import BlastHitFiltering
        hit_filtering = BlastHitFiltering(camel)
        SnakemakeUtils.add_pickle_inputs(hit_filtering, input)
        step = Step(rule, hit_filtering, camel, params.running_dir, config)

        # Update parameters
        hit_filtering.update_parameters(output_filename=os.path.join(params.running_dir, params.output_filename))
        if 'blast_filtering_options' in params.db_config:
            hit_filtering.update_parameters(**params.db_config['blast_filtering_options'])
        if 'extra_column' in params.db_config:
           hit_filtering.update_parameters(
               extra_column_name=params.db_config['extra_column']['name'],
               extra_column_key=params.db_config['extra_column']['key'])

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
        from camel.app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter(camel)
        SnakemakeUtils.add_pickle_inputs(blast_formatter, input)
        step = Step(rule, blast_formatter, camel, params.running_dir, config)
        blast_formatter.update_parameters(output_format='0', num_alignments=1000)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blast_formatter, output)

rule Gene_detection_text_alignment_extraction:
    """
    Extracts a text alignment for the selected hits and attaches them to the hit objects.
    """
    input:
        TXT=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'alignment_generation', 'txt.io'),
        VAL_Hits=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_filtering', 'blast-hits.io'),
        INFORMS_db_info=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'informs.io')
    output:
        VAL_Hits=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'alignment_extraction', 'blast-hits.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'alignment_extraction')
    run:
        from camel.app.tools.pipelines.genedetection.alignmentextractor import AlignmentExtractor
        alignment_extractor = AlignmentExtractor(camel)
        SnakemakeUtils.add_pickle_inputs(alignment_extractor, input)
        step = Step(rule, alignment_extractor, camel, params.running_dir, config)
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
    Input is a pickled dictionary with ToolIO files with either 'FASTQ_PE' or 'FASTQ_SE' as key.
    If paired end input is provided, the read status ('_1', '_1P') is determined based on the read name. 
    """
    input:
        FASTQ = os.path.join(config['working_dir'], INPUT_GENE_DETECTION_FASTQ),
        FASTA = os.path.join(config['working_dir'], 'gene_detection', '{db}', 'db_manager', 'fasta-clust.io')
    output:
        TSV = os.path.join(config['working_dir'], 'gene_detection', '{db}', 'srst2', 'tsv.io'),
        INFORMS = os.path.join(config['working_dir'], 'gene_detection', '{db}', 'srst2', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'srst2'),
        db_config = lambda wildcards: config['gene_detection'][wildcards.db]
    threads: 4
    run:
        from camel.app.tools.srst2.srst2gene import Srst2Gene
        srst2 = Srst2Gene(camel)
        input_files = SnakemakeUtils.load_object(input.FASTQ)
        SnakemakeUtils.add_pickle_input(srst2, 'FASTA', input.FASTA)
        srst2.add_input_files({'FASTQ_PE' if len(input_files) == 2 else 'FASTQ_SE': input_files})
        step = Step(rule, srst2, camel, params.running_dir, config)

        # Update parameters
        srst2.update_parameters(threads=threads)
        if len(input_files) == 2:
            fwd_read_path = input_files[0].path
            fwd_designator, rev_designator = SequenceTypingUtils.determine_read_status(fwd_read_path)
            srst2.update_parameters(forward_designator=fwd_designator, reverse_designator=rev_designator)
        if 'srst2_options' in params.db_config:
            srst2.update_parameters(**params.db_config['srst2_options'])

        # Run tool
        step.run_step()
        SnakemakeUtils.dump_object(srst2.informs, output.INFORMS)
        if 'TSV' in srst2.tool_outputs:
            SnakemakeUtils.dump_tool_output(srst2, 'TSV', output.TSV)
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
        output_filename = lambda wildcards: 'hits-{}-{}.tsv'.format(
            FileSystemHelper.make_valid(str(config['sample_name'])),
            FileSystemHelper.make_valid(wildcards.db)),
        db_config = lambda wildcards: config['gene_detection'][wildcards.db]
    run:
        from camel.app.tools.pipelines.genedetection.srst2hitextractor import SRST2HitExtractor
        extractor = SRST2HitExtractor(camel)
        step = Step(rule, extractor, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(extractor, input)
        extractor.update_parameters(output_filename=os.path.join(params.running_dir, params.output_filename))

        # Add column with additional metadata
        if 'extra_column' in params.db_config:
            extractor.update_parameters(
                extra_column_name=params.db_config['extra_column']['name'],
                extra_column_key=params.db_config['extra_column']['key'])

        step.run_step()
        SnakemakeUtils.dump_tool_outputs(extractor, output)

rule Gene_detection_get_hits:
    """
    Retrieves the hits from the blastn / SRST2 detection method based on the config
    """
    input:
        hits_blast=lambda wildcards: os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_HITS_BLAST.format(db=wildcards.db)) if GeneDetectionUtils.get_detection_method_key(config, wildcards.db) == 'blast' else [],
        hits_srst2=lambda wildcards: os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_HITS_SRST2.format(db=wildcards.db)) if GeneDetectionUtils.get_detection_method_key(config, wildcards.db) == 'srst2' else [],
        tsv_blast=lambda wildcards: os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_TSV_BLAST.format(db=wildcards.db)) if GeneDetectionUtils.get_detection_method_key(config, wildcards.db) == 'blast' else [],
        tsv_srst2=lambda wildcards: os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_TSV_SRST2.format(db=wildcards.db)) if GeneDetectionUtils.get_detection_method_key(config, wildcards.db) == 'srst2' else [],
        informs_blast=lambda wildcards: os.path.join(config['working_dir'], 'gene_detection', '{db}', 'blastn', 'informs.io') if GeneDetectionUtils.get_detection_method_key(config, wildcards.db) == 'blast' else [],
        informs_srst2=lambda wildcards: os.path.join(config['working_dir'], 'gene_detection', '{db}', 'srst2', 'informs.io') if GeneDetectionUtils.get_detection_method_key(config, wildcards.db) == 'srst2' else []
    output:
        VAL_Hits = os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_ALL_HITS),
        TSV = os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_selection', 'selected-tsv.io'),
        INFORMS = os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_INFORMS)
    run:
        shutil.copyfile(input.hits_blast if len(input.hits_blast) > 0 else input.hits_srst2, output.VAL_Hits)
        shutil.copyfile(input.tsv_blast if len(input.tsv_blast) > 0 else input.tsv_srst2, output.TSV)
        shutil.copyfile(input.informs_blast if len(input.informs_blast) > 0 else input.informs_srst2, output.INFORMS)

rule Gene_detection_get_column_names:
    output:
        INFORMS_columns=os.path.join(config['working_dir'], OUTPUT_GENE_DETECTION_COLUMNS)
    params:
        detection_method=lambda wildcards: GeneDetectionUtils.get_detection_method_key(config, wildcards.db),
        extra_column=lambda wildcards: config['gene_detection'][wildcards.db].get('extra_column')
    run:
        from camel.app.components.genedetection.genedetectionblasthit import GeneDetectionBlastHit
        from camel.app.components.genedetection.genedetectionsrst2hit import GeneDetectionSRST2Hit
        if params.detection_method == 'blast':
            columns = GeneDetectionBlastHit.get_column_names_html(params.extra_column[0] if params.extra_column is not None else None)
        elif params.detection_method == 'srst2':
            columns = GeneDetectionSRST2Hit.get_column_names_html(params.extra_column[0] if params.extra_column is not None else None)
        else:
            raise ValueError(f"Invalid detection method: {params.detection_method}")
        SnakemakeUtils.dump_object(columns, output.INFORMS_columns)

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
        sample_name=config['sample_name'],
        config_data=lambda wildcards: config['gene_detection'][wildcards.db],
    run:
        from camel.app.tools.pipelines.genedetection.htmlreportergenedetection import HtmlReporterGeneDetection
        reporter = HtmlReporterGeneDetection(camel)
        step = Step(rule, reporter, camel, params.running_dir, config)
        if 'force_detection_method' in params.config_data:
            reporter.update_parameters(forced_detection_method = params.config_data['force_detection_method'])
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
        INFORMS_hits=os.path.join(config['working_dir'], 'gene_detection', '{db}', 'hit_selection', 'selected-hits.io'),
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
