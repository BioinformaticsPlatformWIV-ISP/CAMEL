from app.components.html.htmlelement import HtmlElement
from app.io.tooliovalue import ToolIOValue


def get_loci(folder):
    """
    Returns the loci from the given database.
    :param folder: Database folder.
    :return: Loci from database
    """
    # for testing, only first 10 genes
    # return [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d)) and not d.startswith('.')][0:10]
    # for real, all genes
    return [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d)) and not d.startswith('.')]


def has_profiles(scheme):
    """
    Returns True if the scheme has profiles.
    :param scheme: Scheme
    :return: True if the scheme has profiles
    """
    folder = SCHEMES.get(scheme)
    return os.path.isfile(os.path.join(folder, 'profiles.tsv'))


SCHEMES = {
    'cgMLST': '/data/sequence_typing/listeria/cgmlst',
    'MLST-Pasteur': '/data/sequence_typing/listeria/mlst/',
    'species_confirmation': '/data/sequence_typing/listeria/species_confirmation/',
    'serogroup': '/data/sequence_typing/listeria/serogroup',
    'virulence': '/data/sequence_typing/listeria/virulence',
    'antibiotic_resistance': '/data/sequence_typing/listeria/antibiotic_resistance',
    'metal_detergent_resistance': '/data/sequence_typing/listeria/metal_detergent_resistance',
}
LOCI = {k: get_loci(SCHEMES.get(k)) for k in SCHEMES.keys()}


rule extract_locus_set_info:
    """
    Extracts the scheme metadata.
    """
    input:
        lambda wildcards: SCHEMES.get(wildcards.scheme)
    output:
        INFORMS = os.path.join(__WORKING_DIR, 'sequence_typing', 'informs-{scheme}.io')
    run:
        from app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
        locus_set_manager = LocusSetManager(camel)
        locus_set_manager.add_input_files({'DIR': [ToolIODirectory(input[0])]})
        locus_set_manager.run()
        SnakemakeUtils.dump_object(locus_set_manager.informs, output.INFORMS)

rule extract_locus_info:
    """
    Extracts the metadata for a single locus.
    """
    input:
        lambda wildcards: os.path.join(SCHEMES.get(wildcards.scheme), wildcards.locus)
    output:
        FASTA = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'fasta-{locus}.io'),
        INFORMS = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'informs-{locus}.io')
    run:
        from app.tools.pipelines.sequence_typing.locusmanager import LocusManager
        locus_manager = LocusManager(camel)
        locus_manager.add_input_files({'DIR': [ToolIODirectory(input[0])]})
        locus_manager.run('.')
        SnakemakeUtils.dump_tool_outputs(locus_manager, output)

rule pickle_dump_sequence_type_definitions:
    """
    Retrieves the sequence type definitions and converts them to CAMEL IO pickles.
    """
    input:
        lambda wildcards: SCHEMES.get(wildcards.scheme)
    output:
        os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'tsv-profiles.io')
    run:
        SnakemakeUtils.dump_object([ToolIOFile(os.path.join(input[0], 'profiles.tsv'))], output[0])

rule blastn_allele_detection:
    """
    Allele detection using blastn.
    """
    input:
        FASTA = os.path.join(__WORKING_DIR, 'assembly', 'fasta.io'),
        DB_BLAST = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'fasta-{locus}.io'),
        INFORMS_locus = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'informs-{locus}.io')
    output:
        VAL_Hit = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', '{locus}', 'hit-blast.io')
    params:
        working_dir = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', '{locus}')
    threads:
        1
    run:
        from app.tools.blast.blastformatter import BlastFormatter
        from app.tools.blast.blastn import Blastn
        from app.tools.pipelines.sequence_typing.besthitselector import BestHitSelector
        from app.tools.pipelines.sequence_typing.alignmentextractor import AlignmentExtractor

        # Blastn alignment
        blastn = Blastn(camel)
        SnakemakeUtils.add_pickle_input(blastn, 'FASTA', input.FASTA)
        SnakemakeUtils.add_pickle_input(blastn, 'DB_BLAST', input.DB_BLAST)
        blastn.update_parameters(threads=threads)
        blastn.run(params.working_dir)

        # TSV generation
        formatter_tsv = BlastFormatter(camel)
        formatter_tsv.update_parameters(output_format='"7 pident sseqid sseq slen qseqid qstart qend"')
        formatter_tsv.add_input_files({'ASN': blastn.tool_outputs['ASN']})
        formatter_tsv.run(params.working_dir)

        # Best hit selection
        hit_selector = BestHitSelector(camel)
        hit_selector.add_input_files({'TSV': formatter_tsv.tool_outputs['TSV']})
        hit_selector.add_input_informs({'locus': SnakemakeUtils.load_object(input.INFORMS_locus)})
        hit_selector.run(params.working_dir)

        # Text alignment generation
        formatter_text = BlastFormatter(camel)
        formatter_text.update_parameters(output_format='0', num_alignments=1000)
        formatter_text.add_input_files({'ASN': blastn.tool_outputs['ASN']})
        formatter_text.run(params.working_dir)

        # Alignment extraction
        extractor = AlignmentExtractor(camel)
        extractor.add_input_files({'TXT': formatter_text.tool_outputs['TXT'],
                                   'VAL_Hits': hit_selector.tool_outputs['VAL_Hit']})
        extractor.run(params.working_dir)

        # Add the alignment to the hit object
        if len(extractor.tool_outputs['TXT']) > 0:
            best_hit = hit_selector.tool_outputs['VAL_Hit'][0].value
            best_hit.alignment_path = extractor.tool_outputs['TXT'][0].path
        SnakemakeUtils.dump_object(hit_selector.tool_outputs['VAL_Hit'], output.VAL_Hit)

rule srst2_allele_detection:
    """
    Allele detection using SRST2.
    """
    input:
        FASTQ_PE = os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-pe.io'),
        FASTA = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'fasta-{locus}.io'),
        INFORMS_locus = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'informs-{locus}.io')
    output:
        VAL_Hit = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', '{locus}', 'hit-srst2.io')
    params:
        working_dir = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', '{locus}')
    threads:
        4
    run:
        # TODO: Add unpaired reads
        from app.tools.srst2.srst2alleledetector import SRST2AlleleDetector
        detector = SRST2AlleleDetector(camel)
        SnakemakeUtils.add_pickle_inputs(detector, input)
        detector.update_parameters(threads=threads, forward_designator='1P', reverse_designator='2P')
        detector.run(params.working_dir)
        SnakemakeUtils.dump_tool_outputs(detector, output)

rule combine_hits:
    """
    Combines hits from one scheme into a list. This also ensures that the allele detection is done using blast or
    using SRST2 depending on the config.
    """
    input:
        hits_blast = lambda wildcards: expand(os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', '{locus}', 'hit-blast.io'),
                                              locus=LOCI.get(wildcards.scheme), scheme=wildcards.scheme) if config['detection_method'] == 'fast' else [],
        hits_srst2 = lambda wildcards: expand(os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', '{locus}', 'hit-srst2.io'),
                                              locus=LOCI.get(wildcards.scheme), scheme=wildcards.scheme) if config['detection_method'] == 'normal' else []
    output:
        os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'all-hits.io')
    run:
        list_of_hits = []
        for pickle in input.hits_blast + input.hits_srst2:
            list_of_hits += SnakemakeUtils.load_object(pickle)
        SnakemakeUtils.dump_object(list_of_hits, output[0])

rule combine_loci:
    """
    Combines the output of all loci of a scheme into a tabular output file.
    """
    input:
        hits = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'all-hits.io'),
        INFORMS_Loci = lambda wildcards: expand(os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'informs-{locus}.io'),
                                                locus=LOCI.get(wildcards.scheme), scheme=wildcards.scheme)
    output:
        TSV = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'tsv-combined.io')
    params:
        working_dir = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}'),
        output_filename = lambda wildcards: os.path.join(__WORKING_DIR, 'sequence_typing', 'typing-{}-{}.tsv'.format(
            FileSystemHelper.make_valid(config['sample_name']),
            FileSystemHelper.make_valid(wildcards.scheme)))
    run:
        from app.tools.pipelines.sequence_typing.allelecombiner import AlleleCombiner
        list_of_informs = []
        for pickle in input.INFORMS_Loci:
            list_of_informs.append(SnakemakeUtils.load_object(pickle))
        combiner = AlleleCombiner(camel)
        combiner.update_parameters(output_filename=params.output_filename)
        combiner.add_input_files({'VAL_Hits': SnakemakeUtils.load_object(input.hits)})
        combiner.add_input_informs({'Loci': list_of_informs})
        combiner.run(params.working_dir)
        SnakemakeUtils.dump_tool_outputs(combiner, output)

rule detect_sequence_type:
    """
    Detects the sequence type based on the detected alleles.
    """
    input:
        VAL_Hits = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'all-hits.io'),
        TSV = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'tsv-profiles.io')
    output:
        INFORMS = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'sequence_type-informs.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}')
    run:
        from app.tools.pipelines.sequence_typing.listeria.sequencetypedetector import SequenceTypeDetector
        sequence_type_detector = SequenceTypeDetector(camel)
        #  SnakemakeUtils.run_tool(sequence_type_detector, input, output, '.')
        SnakemakeUtils.add_pickle_inputs(sequence_type_detector, input)
        step = SnakeStep(rule, sequence_type_detector, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sequence_type_detector, output)


rule create_report:
    """
    Creates a report with the sequence typing output.
    """
    input:
        TSV = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'tsv-combined.io'),
        INFORMS_Scheme = os.path.join(__WORKING_DIR, 'sequence_typing', 'informs-{scheme}.io'),
        INFORMS_ST = lambda wildcards: os.path.join(__WORKING_DIR, 'sequence_typing', '{}', 'sequence_type-informs.io').format(
            wildcards.scheme) if has_profiles(wildcards.scheme) else [],
        VAL_Hits = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'all-hits.io'),
    output:
        VAL_HTML = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'report_sequence_typing', 'html.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'report_sequence_typing'),
    run:
        from app.tools.pipelines.sequence_typing.htmlreportertyping import HtmlReporterTyping
        reporter = HtmlReporterTyping(camel)
        if len(input.INFORMS_ST) != 0:
            reporter.add_input_informs({'ST': SnakemakeUtils.load_object(input.INFORMS_ST)})
        SnakemakeUtils.add_pickle_inputs(reporter, input, ['TSV', 'INFORMS_Scheme', 'VAL_Hits'])
        reporter.run(params.running_dir)
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule combine_sequence_typing_reports:
    """
    Combines the reports of all sequence typing schemes.
    """
    input:
        VAL_HTML = expand(os.path.join(__WORKING_DIR, 'sequence_typing', '{scheme}', 'report_sequence_typing', 'html.io'),
                          scheme=config.get('sequence_typing', None))
    output:
        VAL_HTML = os.path.join(__WORKING_DIR, 'report_sequence_typing', 'html.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'report_sequence_typing'),
        output_dir = config['output_dir']
    run:
        typing_module = HtmlElement('div')
        if len(input.VAL_HTML) > 0:
            typing_module.add_module_header('Sequence Typing')
            for pickle in input.VAL_HTML:
                report_section = SnakemakeUtils.load_object(pickle)[0].value
                report_section.copy_files(params.output_dir)
                typing_module.add_html_object(report_section)
        SnakemakeUtils.dump_object([ToolIOValue(typing_module)], output[0])
