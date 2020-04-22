from pathlib import Path


from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.genedetection.genedetectionutils import GeneDetectionUtils
from camel.app.components.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import gene_detection

camel = Camel.get_instance()


rule gene_detection_db_manager:
    """
    Retrieves the FASTA file and the metadata from a database folder.
    """
    output:
        FASTA = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'db_manager' / 'fasta.io',
        FASTA_clustered = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'db_manager' / 'fasta-clust.io',
        INFORMS = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'db_manager' / 'informs.io'
    params:
        db_path = lambda wildcards: config['gene_detection'][wildcards.db]['path'],
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db
    run:
        from camel.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.pipelines.genedetection.dbmanager import DBManager
        db_manager = DBManager(camel)
        db_manager.add_input_files({'DIR': [ToolIODirectory(params.db_path)]})
        step = Step(rule, db_manager, camel, params.running_dir, config, wildcards)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(db_manager, output)

##########
# BLASTN #
##########
rule gene_detection_blastn:
    """
    Performs local alignment using Blastn+.
    """
    input:
        FASTA = Path(config['working_dir']) / gene_detection.INPUT_GENE_DETECTION_FASTA,
        DB_BLAST = rules.gene_detection_db_manager.output.FASTA_clustered
    output:
        ASN = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'blastn' / 'asn.io',
        INFORMS = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'blastn' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'blastn',
        task = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('blastn', {}).get('task', 'megablast')
    run:
        from camel.app.tools.blast.blastn import Blastn
        blastn = Blastn(camel)
        SnakemakeUtils.add_pickle_inputs(blastn, input)
        step = Step(rule, blastn, camel, params.running_dir, config, wildcards)
        blastn.update_parameters(threads=1, task=params.task)
        step.run_step()
        blastn.informs['Task'] = params.task
        SnakemakeUtils.dump_tool_outputs(blastn, output)

rule gene_detection_tsv_generation:
    """
    Generates tabular output format to extract hit statistics.
    """
    input:
        ASN = rules.gene_detection_blastn.output.ASN
    output:
        TSV = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'tsv_generation' / 'tsv.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'tsv_generation'
    run:
        from camel.app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter(camel)
        SnakemakeUtils.add_pickle_inputs(blast_formatter, input)
        step = Step(rule, blast_formatter, camel, params.running_dir, config, wildcards)
        blast_formatter.update_parameters(output_format='"7 pident sseqid sseq slen qseqid qstart qend score"')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blast_formatter, output)

rule gene_detection_hit_filtering:
    """
    Filters hits based on percent identity and query coverage.
    Extracts the hit information based on the database metadata.
    """
    input:
        TSV = rules.gene_detection_tsv_generation.output.TSV,
        INFORMS_db_info = rules.gene_detection_db_manager.output.INFORMS,
        INFORMS_blastn = rules.gene_detection_blastn.output.INFORMS
    output:
        VAL_Hits = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'hit_filtering' / 'blast-hits.io',
        INFORMS = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS_METHOD).format(db='{db}', method='blast'),
        TSV = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_TABULAR_METHOD).format(db='{db}', method='blast')
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'hit_filtering',
        output_filename = lambda wildcards: 'hits-{}-{}.tsv'.format(
            FileSystemHelper.make_valid(config['sample_name']), FileSystemHelper.make_valid(wildcards.db)),
        db_config = lambda wildcards: config['gene_detection'][wildcards.db],
        min_percent_identity = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('blast', {}).get('min_percent_identity', 90),
        min_coverage = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('blast', {}).get('min_coverage', 60),
    run:
        from camel.app.tools.pipelines.genedetection.blasthitfiltering import BlastHitFiltering
        hit_filtering = BlastHitFiltering(camel)
        SnakemakeUtils.add_pickle_inputs(hit_filtering, input)
        step = Step(rule, hit_filtering, camel, params.running_dir, config, wildcards)

        # Update parameters
        hit_filtering.update_parameters(
            output_filename=str(Path(params.running_dir) / params.output_filename),
            min_percent_identity=str(params.min_percent_identity),
            min_coverage=str(params.min_coverage)
        )
        if params.db_config.get('metadata') is not None:
            hit_filtering.update_parameters(
                extra_column_name=params.db_config['metadata']['name'],
                extra_column_key=params.db_config['metadata']['key'])

        # Run tool
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(hit_filtering, output)

        # Add the informs from the filtering to the existing ones with the blastn command
        informs = SnakemakeUtils.load_object(input.INFORMS_blastn)
        for key, value in hit_filtering.informs.items():
            if key.startswith('_'):
                continue
            informs[key] = value
        SnakemakeUtils.dump_object(informs, output.INFORMS)

rule gene_detection_text_alignment_generation:
    """
    Generates alignments in the text format.
    """
    input:
        ASN = rules.gene_detection_blastn.output.ASN
    output:
        TXT = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'alignment_generation' / 'txt.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'alignment_generation'
    run:
        from camel.app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter(camel)
        SnakemakeUtils.add_pickle_inputs(blast_formatter, input)
        step = Step(rule, blast_formatter, camel, params.running_dir, config, wildcards)
        blast_formatter.update_parameters(output_format='0', num_alignments=1000)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blast_formatter, output)

rule gene_detection_text_alignment_extraction:
    """
    Extracts a text alignment for the selected hits and attaches them to the hit objects.
    """
    input:
        TXT = rules.gene_detection_text_alignment_generation.output.TXT,
        VAL_Hits = rules.gene_detection_hit_filtering.output.VAL_Hits
    output:
        VAL_Hits = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_HITS_METHOD).format(method='blast', db='{db}')
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'alignment_generation'
    run:
        from camel.app.tools.pipelines.genedetection.alignmentextractor import AlignmentExtractor
        alignment_extractor = AlignmentExtractor(camel)
        SnakemakeUtils.add_pickle_inputs(alignment_extractor, input)
        step = Step(rule, alignment_extractor, camel, params.running_dir, config, wildcards)
        step.run_step()
        hits_with_alignment = []
        for io_value, alignment in zip(
                SnakemakeUtils.load_object(input.VAL_Hits), alignment_extractor.tool_outputs['TXT']):
            io_value.value.alignment_path = Path(alignment.path)
            hits_with_alignment.append(io_value)
        SnakemakeUtils.dump_object(hits_with_alignment, output.VAL_Hits)

#########
# SRST2 #
#########
rule gene_detection_srst2:
    """
    Read-mapping based gene detection using SRST2.
    Input is a pickled dictionary with ToolIO files with either 'FASTQ_PE' or 'FASTQ_SE' as key.
    If paired end input is provided, the read status ('_1', '_1P') is determined based on the read name. 
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
        FASTA = rules.gene_detection_db_manager.output.FASTA_clustered
    output:
        TSV = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'srst2' / 'tsv-srst2.io',
        INFORMS = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS_METHOD).format(db='{db}', method='srst2')
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'srst2',
        db_config = lambda wildcards: config['gene_detection'][wildcards.db],
        max_divergence = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('srst2', {}).get('max_divergence', 10),
        min_coverage = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('srst2', {}).get('min_coverage', 60),
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE'
    threads: 4
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.srst2.srst2gene import Srst2Gene
        if not params.running_dir.exists():
            params.running_dir.mkdir(parents=True)
        srst2 = Srst2Gene(camel)
        SnakemakeUtils.add_pickle_input(srst2, 'FASTA', input.FASTA)
        fq_input_dict = SnakePipelineUtils.extracts_fq_input(
            input.IO, key_pe='FASTQ_PE', key_se='FASTQ_SE', read_type=params.read_type)
        srst2.add_input_files(fq_input_dict)
        step = Step(rule, srst2, camel, params.running_dir, config, wildcards)

        # Update parameters
        srst2.update_parameters(threads=threads)
        if 'FASTQ_PE' in fq_input_dict:
            fwd_read_path = fq_input_dict['FASTQ_PE'][0].path
            fwd_designator, rev_designator = SequenceTypingUtils.determine_read_status(fwd_read_path)
            srst2.update_parameters(forward_designator=fwd_designator, reverse_designator=rev_designator)
        srst2.update_parameters(max_divergence=str(params.max_divergence), min_coverage=str(params.min_coverage))

        # Run tool
        step.run_step()
        SnakemakeUtils.dump_object(srst2.informs, output.INFORMS)
        if 'TSV' in srst2.tool_outputs:
            SnakemakeUtils.dump_tool_output(srst2, 'TSV', output.TSV)
        else:
            SnakemakeUtils.dump_object([], output.TSV)

rule gene_detection_srst2_hit_extraction:
    """
    Extracts hits from the SRST2 output.
    """
    input:
        TSV = rules.gene_detection_srst2.output.TSV,
        INFORMS_db = rules.gene_detection_db_manager.output.INFORMS
    output:
        VAL_Hits = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_HITS_METHOD).format(db='{db}', method='srst2'),
        TSV = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_TABULAR_METHOD).format(db='{db}', method='srst2')
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'srst2',
        output_filename = lambda wildcards: 'hits-{}-{}.tsv'.format(
            FileSystemHelper.make_valid(str(config['sample_name'])),
            FileSystemHelper.make_valid(wildcards.db)),
        db_config = lambda wildcards: config['gene_detection'][wildcards.db]
    run:
        from camel.app.tools.pipelines.genedetection.srst2hitextractor import SRST2HitExtractor
        extractor = SRST2HitExtractor(camel)
        step = Step(rule, extractor, camel, params.running_dir, config, wildcards)
        SnakemakeUtils.add_pickle_inputs(extractor, input)
        extractor.update_parameters(output_filename=os.path.join(params.running_dir, params.output_filename))

        # Add column with additional metadata
        if 'metadata' in params.db_config:
            extractor.update_parameters(
                extra_column_name=params.db_config['metadata']['name'],
                extra_column_key=params.db_config['metadata']['key'])

        step.run_step()
        SnakemakeUtils.dump_tool_outputs(extractor, output)

#######
# KMA #
#######
rule gene_detection_kma_get_db:
    """
    Retrieves the database for running KMA.
    """
    input:
        FASTA = rules.gene_detection_db_manager.output.FASTA
    output:
        DB = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'kma' / 'db.io'
    run:
        import json
        from camel.app.io.tooliovalue import ToolIOValue
        fasta_path = Path(SnakemakeUtils.load_object(input.FASTA)[0].path)
        with open(fasta_path.parent / 'db_metadata.txt') as handle:
            metadata = json.load(handle)
        dir_kma = fasta_path.parent / 'kma'
        if not dir_kma.exists():
            raise FileNotFoundError(f"KMA database not found: {dir_kma}")
        kma_path = fasta_path.parent / 'kma' / metadata['name'].lower()
        SnakemakeUtils.dump_object([ToolIOValue(kma_path)], output.DB)

rule gene_detection_kma:
    """
    Runs KMA on a database with the gene detection.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
        DB = rules.gene_detection_kma_get_db.output.DB
    output:
        TSV = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'kma' / 'tsv-kma.io',
        INFORMS = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS_METHOD).format(method='kma', db='{db}')
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'kma'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.kma.kma import KMA
        kma = KMA(camel)
        SnakemakeUtils.add_pickle_input(kma, 'DB', input.DB)
        kma.add_input_files(SnakePipelineUtils.extracts_fq_input(input.IO))
        step = Step(rule, kma, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(kma, output)

rule gene_detection_KMA_hit_extraction:
    """
    Extracts and filters the hits detected by KMA.
    """
    input:
        TSV = rules.gene_detection_kma.output.TSV,
        INFORMS_db = rules.gene_detection_db_manager.output.INFORMS
    output:
        VAL_hits = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_HITS_METHOD).format(method='kma', db='{db}'),
        TSV = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_TABULAR_METHOD).format(method='kma', db='{db}')
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'kma'
    run:
        from camel.app.tools.kma.kmagenedetectionhitextractor import KMAGeneDetectionHitExtractor
        extractor = KMAGeneDetectionHitExtractor(camel)
        SnakemakeUtils.add_pickle_inputs(extractor, input)
        step = Step(rule, extractor, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(extractor, output)

rule gene_detection_get_hits:
    """
    Retrieves the hits from the blastn / SRST2 detection method based on the config
    """
    input:
        informs_db = rules.gene_detection_db_manager.output.INFORMS,
        VAL_hits = lambda wildcards: str(Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_HITS_METHOD).format(db=wildcards.db, method=GeneDetectionUtils.get_detection_method_key(config, wildcards.db))),
        TSV = lambda wildcards: str(Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_TABULAR_METHOD).format(db=wildcards.db, method=GeneDetectionUtils.get_detection_method_key(config, wildcards.db))),
        INFORMS = lambda wildcards: str(Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS_METHOD).format(db=wildcards.db, method=GeneDetectionUtils.get_detection_method_key(config, wildcards.db)))
    output:
        TSV = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'hit_selection' / 'selected-tsv.io',
        VAL_hits = Path(config['working_dir']) / gene_detection.OUTPUT_GENE_DETECTION_ALL_HITS,
        INFORMS = Path(config['working_dir']) / gene_detection.OUTPUT_GENE_DETECTION_INFORMS
    run:
        import shutil
        shutil.copyfile(input.VAL_hits, output.VAL_hits)
        shutil.copyfile(input.TSV, output.TSV)
        # Add a tag for the database to distinguish commands in the output
        informs = SnakemakeUtils.load_object(input.INFORMS)
        informs['_tag'] = SnakemakeUtils.load_object(input.informs_db)['title']
        SnakemakeUtils.dump_object(informs, output.INFORMS)

rule gene_detection_get_column_names:
    """
    Retrieves the column names for the output.
    """
    output:
        INFORMS_columns=Path(config['working_dir']) / gene_detection.OUTPUT_GENE_DETECTION_COLUMNS
    params:
        detection_method = lambda wildcards: GeneDetectionUtils.get_detection_method_key(config, wildcards.db),
        db_config = lambda wildcards: config['gene_detection'][wildcards.db]
    run:
        from camel.app.components.blast.blasthitstatistics import BlastHitStatistics
        from camel.app.components.genedetection.genedetectionblasthit import GeneDetectionBlastHit
        from camel.app.components.genedetection.genedetectionsrst2hit import GeneDetectionSRST2Hit

        # Create empty hit
        if params.detection_method == 'blast':
            empty_hit = GeneDetectionBlastHit('Locus', None, BlastHitStatistics('subject', 0, '', 'query', 0, 0, 0.0))
        elif params.detection_method == 'srst2':
            empty_hit = GeneDetectionSRST2Hit('DB_cluster', 'Locus', None, 0, '', '', 0.0, 0.0, 0.0)
        else:
            raise ValueError(f"Invalid detection method: {params.detection_method}")

        # Add metadata columns
        if 'metadata' in params.db_config:
            empty_hit.add_metadata(params.db_config['metadata']['name'], '')

        # Save column names
        columns = empty_hit.html_column_names
        SnakemakeUtils.dump_object(columns, output.INFORMS_columns)

rule gene_detection_report:
    """
    Creates HTML reports for the gene detection.
    """
    input:
        VAL_Hits = rules.gene_detection_get_hits.output.VAL_hits,
        INFORMS_db_info = rules.gene_detection_db_manager.output.INFORMS,
        INFORMS_detection = rules.gene_detection_get_hits.output.INFORMS,
        TSV = rules.gene_detection_get_hits.output.TSV
    output:
        VAL_HTML = Path(config['working_dir']) / gene_detection.OUTPUT_GENE_DETECTION_REPORT,
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'report',
        sample_name = config['sample_name'],
        config_data = lambda wildcards: config['gene_detection'][wildcards.db],
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.genedetection.htmlreportergenedetection import HtmlReporterGeneDetection
        reporter = HtmlReporterGeneDetection(camel)
        step = Step(rule, reporter, camel, params.running_dir, config, wildcards)
        if 'force_detection_method' in params.config_data:
            reporter.update_parameters(forced_detection_method = params.config_data['force_detection_method'])
        reporter.add_input_files({'SAMPLE_NAME': [ToolIOValue(params.sample_name)]})
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule gene_detection_create_empty_report:
    """
    Creates an empty HTML report for the gene detection.
    """
    input:
        INFORMS_db_info = rules.gene_detection_db_manager.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / gene_detection.OUTPUT_GENE_DETECTION_REPORT_EMPTY
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'report',
    run:
        from camel.app.components.html.htmlreportsection import HtmlReportSection
        db_info = SnakemakeUtils.load_object(input[0])
        section = HtmlReportSection(db_info['title'], 3)
        section.add_paragraph('Analysis disabled')
        SnakemakeUtils.dump_object([ToolIOValue(section)], output[0])

rule gene_detection_dump_summary_info:
    """
    Dumps the summary information from the gene detection in tabular format.
    """
    input:
        INFORMS_hits=Path(config['working_dir']) / 'gene_detection' / '{db}' / 'hit_selection' / 'selected-hits.io',
    output:
        TSV = Path(config['working_dir']) / gene_detection.OUTPUT_GENE_DETECTION_SUMMARY
    run:
        import dataclasses
        informs = SnakemakeUtils.load_object(input.INFORMS_hits)
        hit_info = []
        blast_stats = []
        for hit in informs:
            hit_info.append(hit.value.to_table_row())
        with open(output[0], 'w') as handle:
            handle.write('hits_{}\t{}'.format(wildcards.db, json.dumps(hit_info)))
            handle.write('\n')
