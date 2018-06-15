GENE_DETECTION_WORKING_DIR = os.path.join(__WORKING_DIR, 'gene_detection')
GENE_DETECTION_REPORT = os.path.join(GENE_DETECTION_WORKING_DIR, 'report-html.io')
GENE_DETECTION_DB_SUMMARY = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'report', 'summary.tsv')

rule database_manager:
    """
    Retrieves the FASTA file and the metadata from a database folder.
    """
    output:
        FASTA = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'database_manager', 'fasta.io'),
        INFORMS = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'database_manager', 'informs.io')
    params:
        running_dir = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}')
    run:
        from camel.app.tools.pipelines.gene_detection.dbmanager import DBManager
        db_manager = DBManager(camel)
        db_manager.add_input_files({'DIR': [ToolIODirectory(ToolIODb(wildcards.db).path)]})
        SnakemakeUtils.add_pickle_inputs(db_manager, input)
        step = Step(rule, db_manager, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(db_manager, output)

rule blastn:
    """
    Performs local alignment using Blastn+.
    """
    input:
        FASTA = FASTA_ASSEMBLY,
        DB_BLAST = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'database_manager', 'fasta.io')
    output:
        ASN = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'blastn', 'asn.io')
    params:
        running_dir = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'blastn')
    run:
        from camel.app.tools.blast.blastn import Blastn
        blastn = Blastn(camel)
        SnakemakeUtils.add_pickle_inputs(blastn, input)
        step = Step(rule, blastn, camel, params.running_dir, config)
        blastn.update_parameters(threads=1)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blastn, output)

rule tsv_generation:
    """
    Generates tabular output format to extract hit statistics.
    """
    input:
        ASN = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'blastn', 'asn.io')
    output:
        TSV = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'tsv_generation', 'tsv.io')
    params:
        running_dir = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'tsv_generation')
    run:
        from camel.app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter(camel)
        SnakemakeUtils.add_pickle_inputs(blast_formatter, input)
        step = Step(rule, blast_formatter, camel, params.running_dir, config)
        blast_formatter.update_parameters(output_format='"7 pident sseqid sseq slen qseqid qstart qend"')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blast_formatter, output)

rule hit_filtering:
    """
    Filters hits based on percent identity and query coverage.
    Extracts the hit information based on the database metadata.
    """
    input:
        TSV = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'tsv_generation', 'tsv.io'),
        INFORMS_db_info = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'database_manager', 'informs.io')
    output:
        VAL_Hits = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'hit_filtering', 'blast-hits.io'),
        TSV = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'hit_filtering', 'tsv-filtered.io')
    params:
        running_dir = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'hit_filtering'),
        filtering_method=lambda wildcards: config['gene_detection'][wildcards.db].get('filtering_method'),
        min_percent_identity = lambda wildcards: config['gene_detection'][wildcards.db].get('min_percent_identity'),
        min_coverage = lambda wildcards: config['gene_detection'][wildcards.db].get('min_coverage'),
        extra_column = lambda wildcards: config['gene_detection'][wildcards.db].get('extra_column'),
        output_filename = lambda wildcards: os.path.join(GENE_DETECTION_WORKING_DIR, wildcards.db, 'hits-{}-{}.tsv'.format(
            FileSystemHelper.make_valid(config['sample_name']),
            FileSystemHelper.make_valid(wildcards.db)))
    run:
        from camel.app.tools.pipelines.gene_detection.hitfiltering import HitFiltering
        hit_filtering = HitFiltering(camel)
        SnakemakeUtils.add_pickle_inputs(hit_filtering, input)
        step = Step(rule, hit_filtering, camel, params.running_dir, config)

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

rule text_alignment_generation:
    """
    Generates alignments in the text format.
    """
    input:
        ASN = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'blastn', 'asn.io')
    output:
        TXT = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'alignment_generation', 'txt.io')
    params:
        running_dir = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'alignment_generation')
    run:
        from camel.app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter(camel)
        SnakemakeUtils.add_pickle_inputs(blast_formatter, input)
        step = Step(rule, blast_formatter, camel, params.running_dir, config)
        blast_formatter.update_parameters(output_format='0', num_alignments=1000)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blast_formatter, output)

rule text_alignment_extraction:
    """
    Extracts a text alignment for the selected hits and attaches them to the hit objects.
    """
    input:
        TXT = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'alignment_generation', 'txt.io'),
        VAL_Hits = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'hit_filtering', 'blast-hits.io'),
        INFORMS_db_info=os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'database_manager', 'informs.io')
    output:
        VAL_Hits = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'alignment_extraction', 'blast-hits.io')
    params:
        running_dir = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'alignment_extraction')
    run:
        from camel.app.tools.pipelines.gene_detection.alignmentextractor import AlignmentExtractor
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

rule srst2_gene_detection:
    """
    Read-mapping based gene detection using SRST2.
    """
    input:
        FASTQ_PE = TRIMMED_READS_PE,
        FASTA = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'database_manager', 'fasta.io'),
    output:
        TSV = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'srst2', 'tsv.io')
    params:
        running_dir = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'srst2')
    threads:
        8
    run:
        from camel.app.tools.srst2.srst2gene import Srst2Gene
        srst2 = Srst2Gene(camel)
        SnakemakeUtils.add_pickle_inputs(srst2, input)
        step = Step(rule, srst2, camel, params.running_dir, config)
        srst2.update_parameters(threads=threads, forward_designator='1P', reverse_designator='2P')
        step.run_step()
        if 'TSV' in srst2.tool_outputs:
            SnakemakeUtils.dump_tool_outputs(srst2, output)
        else:
            SnakemakeUtils.dump_object([], output.TSV)

rule srst2_hit_extraction:
    """
    Extracts hits from the SRST2 output.
    """
    input:
        TSV = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'srst2', 'tsv.io'),
        INFORMS_db = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'database_manager', 'informs.io'),
    output:
        VAL_Hits = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'srst2', 'srst2-hits.io'),
        TSV = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'srst2', 'tsv-srst2.io')
    params:
        running_dir = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'srst2'),
        output_filename = lambda wildcards: os.path.join(GENE_DETECTION_WORKING_DIR, wildcards.db, 'hits-{}-{}.tsv'.format(
            FileSystemHelper.make_valid(config['sample_name']),
            FileSystemHelper.make_valid(wildcards.db)))
    run:
        from camel.app.tools.pipelines.gene_detection.srst2hitextractor import SRST2HitExtractor
        extractor = SRST2HitExtractor(camel)
        step = Step(rule, extractor, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(extractor, input)
        extractor.update_parameters(output_filename=params.output_filename)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(extractor, output)

rule gene_detection_get_hits:
    """
    Retrieves the hits from the blastn / SRST2 detection method based on the config
    """
    input:
        hits_blast = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'alignment_extraction', 'blast-hits.io') if config[
            'detection_method'] == 'fast' else [],
        hits_srst2 = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'srst2', 'srst2-hits.io') if config[
            'detection_method'] == 'normal' else [],
        tsv_blast = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'hit_filtering', 'tsv-filtered.io') if config[
            'detection_method'] == 'fast' else [],
        tsv_srst2 = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'srst2', 'tsv-srst2.io') if config[
            'detection_method'] == 'normal' else []
    output:
        VAL_Hits = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'selected-hits.io'),
        TSV = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'selected-tsv.io')
    run:
        if len(input.hits_blast) > 0:
            shutil.copyfile(input.hits_blast, output.VAL_Hits)
        else:
            shutil.copyfile(input.hits_srst2, output.VAL_Hits)
        if len(input.tsv_blast) > 0:
            shutil.copyfile(input.tsv_blast, output.TSV)
        else:
            shutil.copyfile(input.tsv_srst2, output.TSV)

rule report_gene_detection:
    """
    Creates HTML reports for the gene detection.
    """
    input:
        VAL_Hits = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'selected-hits.io'),
        INFORMS_db_info = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'database_manager', 'informs.io'),
        TSV = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'selected-tsv.io')
    output:
        VAL_HTML = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'report', 'html.io')
    params:
        running_dir = os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'report'),
        sample_name = config['sample_name']
    run:
        from camel.app.tools.pipelines.gene_detection.htmlreportergenedetection import HtmlReporterGeneDetection
        reporter = HtmlReporterGeneDetection(camel)
        step = Step(rule, reporter, camel, params.running_dir, config)
        reporter.add_input_files({'SAMPLE_NAME': [ToolIOValue(params.sample_name)]})
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule summary_gene_detection:
    """
    Creates a tabular summary for gene detection on a {db}.
    """
    input:
        INFORMS_hits=os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'selected-hits.io'),
    output:
        GENE_DETECTION_DB_SUMMARY
    params:
        running_dir=os.path.join('gene_detection', '{db}', 'report')
    run:
        informs = SnakemakeUtils.load_object(input.INFORMS_hits)
        hit_info = []
        for hit in informs:
            hit_info.append(hit.value.to_table_row().split('\t'))
        with open(output[0], 'w') as handle:
            handle.write('hits_{}\t{}'.format(wildcards.db, json.dumps(hit_info)))
            handle.write('\n')

rule combine_gene_detection_reports:
    """
    Combines the reports from the different databases.
    """
    input:
        HTML_reports = expand(os.path.join(GENE_DETECTION_WORKING_DIR, '{db}', 'report', 'html.io'), db=config['gene_detection']),
    output:
        GENE_DETECTION_REPORT
    params:
        output_dir = config['output_dir'],
        dbs = config['gene_detection']
    run:
        gene_detection_module = HtmlElement('div')
        gene_detection_db_sections = {
            'Resistance Characterization': ['resfinder' ,'card', 'arg_annot'],
            'Virulence Detection': ['virulencefinder_listeria'],
            'Plasmid Replicon Detection': ['plasmidfinder_grampositive']
        }
        section_reports = {key: [] for key in gene_detection_db_sections.keys()}
        for section, dbs in gene_detection_db_sections.items():
            for db in dbs:
                if db in params.dbs:
                    section_reports[section].append(os.path.join(GENE_DETECTION_WORKING_DIR, db, 'report', 'html.io'))

        for section, reports in section_reports.items():
            if len(reports) > 0:
                gene_detection_module.add_module_header(section)
                for pickle in reports:
                    report_section = SnakemakeUtils.load_object(pickle)[0].value
                    report_section.copy_files(params.output_dir)
                    gene_detection_module.add_html_object(report_section)
        SnakemakeUtils.dump_object([ToolIOValue(gene_detection_module)], output[0])
