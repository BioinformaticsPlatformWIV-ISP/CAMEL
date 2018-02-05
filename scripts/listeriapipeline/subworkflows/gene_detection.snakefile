from app.components.html.htmlelement import HtmlElement
from app.components.html.htmlreportsection import HtmlReportSection

class GeneDatabase(object):

    """
    This class represents a gene database.
    """

    def __init__(self, name, path, path_clustered, extra_column, pi_threshold):
        """
        Initializes a gene database.
        :param name: Database name
        :param path: Path to the database directory.
        :param path_clustered: Path to the clustered database directory
        :param extra_column: Extra column containing metadata (name, key)
        :param pi_threshold: Percent identity threshold
        """
        self.name = name
        self.path = path
        self.path_clustered = path_clustered
        self.extra_column = extra_column
        self.pi_threshold = pi_threshold


# Database configuration (possibly move?)
DATABASES = {
    'resistance': [
        GeneDatabase('ARG-ANNOT',
                     '/data/sequence_db/ARG-ANNOT',
                     '/data/srst2/gene_db/ARG-ANNOT-clustered_80',
                     None, 90),
        GeneDatabase('CARD',
                     '/data/sequence_db/CARD',
                     '/data/srst2/gene_db/CARD-clustered_80',
                     None, 90),
        GeneDatabase('ResFinder',
                     '/data/sequence_db/ResFinder',
                     '/data/srst2/gene_db/ResFinder-clustered_80',
                     None, 98)
    ],
    'virulence': [
        GeneDatabase('VirulenceFinder',
                     '/data/sequence_db/VirulenceFinder-Listeria/',
                     '/data/srst2/gene_db/VirulenceFinder-Listeria-clustered_80/',
                     ('Protein function', 'protein_function',), 90),
    ],
    'plasmid': [
        GeneDatabase('Gram_positive',
                     '/data/sequence_db/PlasmidFinder-gram_positive/',
                     '/data/srst2/gene_db/PlasmidFinder-gram_positive-clustered_80/',
                     ['Notes', 'notes'], 95)
    ]
}

def get_db_keys(db_group):
    """
    Returns the database keys for the databases belong to a db group.
    :param db_group: The database group defined in DATABASES (e.g. 'resistance')
    :return: Database keys
    """
    return [d.name for d in DATABASES[db_group] if d.name in config['gene_detection']]

def get_database(name):
    """
    Returns the database with the given name.
    :param name: Database name
    :return: Database
    """
    for db_type in DATABASES:
        for db in DATABASES[db_type]:
            if db.name == name:
                return db
    raise ValueError("Database '{}' not found".format(name))



rule database_manager:
    """
    Retrieves the FASTA file and the metadata from a database folder.
    """
    output:
        FASTA = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'database_manager', 'fasta.io'),
        INFORMS = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'database_manager', 'informs.io')
    params:
        db_path = lambda wildcards: get_database(wildcards.db).path,
        running_dir = os.path.join(__WORKING_DIR, 'gene_detection', '{db}')
    run:
        from app.tools.pipelines.gene_detection.dbmanager import DBManager
        db_manager = DBManager(camel)
        db_manager.add_input_files({'DIR': [ToolIODirectory(params.db_path)]})
        SnakemakeUtils.add_pickle_inputs(db_manager, input)
        step = SnakeStep(rule, db_manager, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(db_manager, output)

rule blastn:
    """
    Performs local alignment using Blastn+.
    """
    input:
        FASTA = os.path.join(__WORKING_DIR, 'assembly', 'fasta.io'),
        DB_BLAST = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'database_manager', 'fasta.io')
    output:
        ASN = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'blastn', 'asn.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'blastn')
    run:
        from app.tools.blast.blastn import Blastn
        blastn = Blastn(camel)
        SnakemakeUtils.add_pickle_inputs(blastn, input)
        step = SnakeStep(rule, blastn, camel, params.running_dir, config)
        blastn.update_parameters(threads=1)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blastn, output)

rule tsv_generation:
    """
    Generates tabular output format to extract hit statistics.
    """
    input:
        ASN = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'blastn', 'asn.io')
    output:
        TSV = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'tsv_generation', 'tsv.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'tsv_generation')
    run:
        from app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter(camel)
        SnakemakeUtils.add_pickle_inputs(blast_formatter, input)
        step = SnakeStep(rule, blast_formatter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blast_formatter, output)

rule hit_filtering:
    """
    Filters hits based on percent identity and query coverage.
    Extracts the hit information based on the database metadata.
    """
    input:
        TSV = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'tsv_generation', 'tsv.io'),
        INFORMS_db_info = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'database_manager', 'informs.io')
    output:
        VAL_Hits = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'hit_filtering', 'blast-hits.io'),
        TSV = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'hit_filtering', 'tsv-filtered.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'hit_filtering'),
        percent_identity = lambda wildcards: get_database(wildcards.db).pi_threshold,
        extra_column = lambda wildcards: get_database(wildcards.db).extra_column,
        output_filename = lambda wildcards: os.path.join(__WORKING_DIR, 'gene_detection', wildcards.db, 'hits-{}-{}.tsv'.format(
            FileSystemHelper.make_valid(config['sample_name']),
            FileSystemHelper.make_valid(wildcards.db)))
    run:
        from app.tools.pipelines.gene_detection.hitfiltering import HitFiltering
        hit_filtering = HitFiltering(camel)
        SnakemakeUtils.add_pickle_inputs(hit_filtering, input)
        step = SnakeStep(rule, hit_filtering, camel, params.running_dir, config)
        hit_filtering.update_parameters(min_percent_identity=params.percent_identity,
                                        output_filename=params.output_filename)
        if params.extra_column is not None:
            hit_filtering.update_parameters(extra_column_name=params.extra_column[0],
                                            extra_column_key=params.extra_column[1])
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(hit_filtering, output)

rule text_alignment_generation:
    """
    Generates alignments in the text format.
    """
    input:
        ASN = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'blastn', 'asn.io')
    output:
        TXT = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'alignment_generation', 'txt.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'alignment_generation')
    run:
        from app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter(camel)
        SnakemakeUtils.add_pickle_inputs(blast_formatter, input)
        step = SnakeStep(rule, blast_formatter, camel, params.running_dir, config)
        blast_formatter.update_parameters(output_format='0', num_alignments=1000)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blast_formatter, output)

rule text_alignment_extraction:
    """
    Extracts a text alignment for the selected hits and attaches them to the hit objects.
    """
    input:
        TXT = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'alignment_generation', 'txt.io'),
        VAL_Hits = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'hit_filtering', 'blast-hits.io')
    output:
        VAL_Hits = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'alignment_extraction', 'blast-hits.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'alignment_extraction')
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

rule get_clustered_db:
    """
    Returns the clustered database that can be used by SRST2.
    """
    output:
        FASTA = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'clustered', 'fasta.io'),
        INFORMS = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'clustered', 'mapping.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'gene_detection', '{db}'),
        db_path = lambda wildcards: get_database(wildcards.db).path_clustered
    run:
        from app.tools.pipelines.gene_detection.dbmanagerclustered import DBManagerClustered
        db_manager = DBManagerClustered(camel)
        step = SnakeStep(rule, db_manager, camel, params.running_dir, config)
        db_manager.add_input_files({'DIR': [ToolIODirectory(params.db_path)]})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(db_manager, output)

rule srst2_gene_detection:
    """
    Read-mapping based gene detection using SRST2.
    """
    input:
        FASTQ_PE = os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-pe.io'),
        FASTA = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'clustered', 'fasta.io'),
    output:
        TSV = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'srst2', 'tsv.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'srst2')
    threads:
        4
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

rule srst2_hit_extraction:
    """
    Extracts hits from the SRST2 output.
    """
    input:
        TSV = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'srst2', 'tsv.io'),
        INFORMS_db_info = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'database_manager', 'informs.io'),
        INFORMS_mapping = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'clustered', 'mapping.io')
    output:
        VAL_Hits = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'srst2', 'srst2-hits.io'),
        TSV = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'srst2', 'tsv-srst2.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'srst2'),
        output_filename = lambda wildcards: os.path.join(__WORKING_DIR, 'gene_detection', wildcards.db, 'hits-{}-{}.tsv'.format(
            FileSystemHelper.make_valid(config['sample_name']),
            FileSystemHelper.make_valid(wildcards.db)))
    run:
        from app.tools.pipelines.gene_detection.srst2hitextractor import SRST2HitExtractor
        extractor = SRST2HitExtractor(camel)
        step = SnakeStep(rule, extractor, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(extractor, input)
        extractor.update_parameters(output_filename=params.output_filename)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(extractor, output)

rule gene_detection_get_hits:
    """
    Retrieves the hits from the blastn / SRST2 detection method based on the config
    """
    input:
        hits_blast = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'alignment_extraction', 'blast-hits.io') if config[
            'detection_method'] == 'fast' else [],
        hits_srst2 = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'srst2', 'srst2-hits.io') if config[
            'detection_method'] == 'normal' else [],
        tsv_blast = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'hit_filtering', 'tsv-filtered.io') if config[
            'detection_method'] == 'fast' else [],
        tsv_srst2 = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'srst2', 'tsv-srst2.io') if config[
            'detection_method'] == 'normal' else []
    output:
        VAL_Hits = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'selected-hits.io'),
        TSV = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'selected-tsv.io')
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
        VAL_Hits = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'selected-hits.io'),
        INFORMS_db_info = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'database_manager', 'informs.io'),
        TSV = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'selected-tsv.io')
    output:
        VAL_HTML = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'report', 'html.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'report'),
        sample_name = config['sample_name']
    run:
        from app.tools.pipelines.gene_detection.htmlreportergenedetection import HtmlReporterGeneDetection
        reporter = HtmlReporterGeneDetection(camel)
        step = SnakeStep(rule, reporter, camel, params.running_dir, config)
        reporter.add_input_files({'SAMPLE_NAME': [ToolIOValue(params.sample_name)]})
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule combine_gene_detection_reports:
    """
    Combines the reports from the different databases.
    """
    input:
        HTML_Res = expand(os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'report', 'html.io'), db=get_db_keys('resistance')),
        HTML_Vir = expand(os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'report', 'html.io'), db=get_db_keys('virulence')),
        HTML_Pla = expand(os.path.join(__WORKING_DIR, 'gene_detection', '{db}', 'report', 'html.io'), db=get_db_keys('plasmid'))
    output:
        os.path.join(__WORKING_DIR, 'report_gene_detection', 'html.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'gene_detection', 'report'),
        output_dir = config['output_dir']
    run:
        gene_detection_module = HtmlElement('div')
        for title, input_html in [
                ('Resistance Characterization', input.HTML_Res),
                ('Virulence Detection', input.HTML_Vir),
                ('Plasmid Replicon Detection', input.HTML_Pla)]:
            if len(input_html) > 0:
                gene_detection_module.add_module_header(title)
                for pickle in input_html:
                    report_section = SnakemakeUtils.load_object(pickle)[0].value
                    report_section.copy_files(params.output_dir)
                    gene_detection_module.add_html_object(report_section)

        SnakemakeUtils.dump_object([ToolIOValue(gene_detection_module)], output[0])
