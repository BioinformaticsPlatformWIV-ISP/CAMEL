"""
This Snakefile can be used to perform sequence typing.
There are two approaches:
- Fast: blastn for nucleotide loci and blastx for protein loci.
- Normal: SRST2 for nucleotide loci and blastx for protein loci.

The sequence typing databases need to be added to the config file under the 'sequence_typing' section.
Each database need to be specified as 'Key': 'Path to DB folder'.
"""
TYPING_WORKING_DIR = os.path.join(__WORKING_DIR, 'sequence_typing')
TYPING_REPORT = os.path.join(TYPING_WORKING_DIR, 'report-html.io')
TYPING_SCHEME_SUMMARY = os.path.join(TYPING_WORKING_DIR, '{scheme}', 'report-summary.io')
SPECIES_CONFIRM_ST_REPORT = os.path.join(TYPING_WORKING_DIR, 'species-confirm-st-report-html.io')

def get_loci(folder):
    """
    Returns the loci from the given database.
    :param folder: Database folder.
    :return: Loci from database
    """
    return sorted([d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d)) and not d.startswith('.')])


def has_profiles(scheme):
    """
    Returns True if the scheme has profiles.
    :param scheme: Scheme
    :return: True if the scheme has profiles
    """
    folder = SCHEMES.get(scheme)
    return os.path.isfile(os.path.join(folder, 'profiles.tsv'))

def get_locus_type(scheme, locus):
    """
    Returns the type of locus ('DNA', 'peptide').
    :param scheme: Scheme (e.g. 'MLST')
    :param locus: Locus (e.g. 'abcZ')
    :return: Locus type
    """
    import json
    locus_directory = os.path.join(SCHEMES[scheme], locus)
    locus_metadata_file = os.path.join(locus_directory, 'locus_metadata.txt')
    if not os.path.isfile(locus_metadata_file):
        raise FileNotFoundError("No metadata found in '{}'".format(locus_directory))
    with open(locus_metadata_file) as handle:
        try:
            return json.load(handle)['type']
        except KeyError:
            raise ValueError("Metadata file does not contain locus type ({})".format(locus_directory))

SCHEMES = config['sequence_typing']
LOCI = {k: get_loci(SCHEMES.get(k)) for k in SCHEMES.keys()}



rule extract_locus_set_info:
    """
    Extracts the scheme metadata.
    """
    input:
        lambda wildcards: SCHEMES.get(wildcards.scheme)
    output:
        INFORMS = os.path.join(TYPING_WORKING_DIR, '{scheme}', 'informs-{scheme}.io')
    params:
        running_dir=os.path.join(TYPING_WORKING_DIR, '{scheme}')
    run:
        from app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
        locus_set_manager = LocusSetManager(camel)
        locus_set_manager.add_input_files({'DIR': [ToolIODirectory(input[0])]})
        step = SnakeStep(rule, locus_set_manager, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_object(locus_set_manager.informs, output.INFORMS)

rule extract_locus_info:
    """
    Extracts the metadata for a single locus.
    """
    input:
        lambda wildcards: os.path.join(SCHEMES.get(wildcards.scheme), wildcards.locus)
    output:
        FASTA = os.path.join(TYPING_WORKING_DIR, '{scheme}', '{locus}', 'fasta-{locus}.io'),
        INFORMS = os.path.join(TYPING_WORKING_DIR, '{scheme}', '{locus}', 'informs-{locus}.io')
    params:
        running_dir=os.path.join(TYPING_WORKING_DIR, '{scheme}')
    run:
        from app.tools.pipelines.sequence_typing.locusmanager import LocusManager
        locus_manager = LocusManager(camel)
        locus_manager.add_input_files({'DIR': [ToolIODirectory(input[0])]})
        step = SnakeStep(rule, locus_manager, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(locus_manager, output)

rule pickle_dump_sequence_type_definitions:
    """
    Retrieves the sequence type definitions and converts them to CAMEL IO pickles.
    """
    input:
        lambda wildcards: SCHEMES.get(wildcards.scheme)
    output:
        os.path.join(TYPING_WORKING_DIR, '{scheme}', 'tsv-profiles.io')
    run:
        SnakemakeUtils.dump_object([ToolIOFile(os.path.join(input[0], 'profiles.tsv'))], output[0])

rule blastn_allele_detection:
    """
    Allele detection using blastn.
    """
    input:
        FASTA = FASTA_ASSEMBLY,
        DB_BLAST = os.path.join(TYPING_WORKING_DIR, '{scheme}', '{locus}', 'fasta-{locus}.io'),
        INFORMS_locus = os.path.join(TYPING_WORKING_DIR, '{scheme}', '{locus}', 'informs-{locus}.io')
    output:
        VAL_Hit = os.path.join(TYPING_WORKING_DIR, '{scheme}', '{locus}', 'hit-blast.io')
    params:
        running_dir = os.path.join(TYPING_WORKING_DIR, '{scheme}', '{locus}')
    threads: 1
    run:
        from app.tools.blast.blastformatter import BlastFormatter
        from app.tools.blast.blastn import Blastn
        from app.tools.blast.blastx import Blastx
        from app.tools.pipelines.sequence_typing.besthitselector import BestHitSelector
        from app.tools.pipelines.sequence_typing.alignmentextractor import AlignmentExtractor

        # Blastn alignment
        if get_locus_type(wildcards.scheme, wildcards.locus) == 'DNA':
            blast = Blastn(camel)
        else:
            blast = Blastx(camel)
        SnakemakeUtils.add_pickle_input(blast, 'FASTA', input.FASTA)
        SnakemakeUtils.add_pickle_input(blast, 'DB_BLAST', input.DB_BLAST)
        blast.update_parameters(threads=threads)
        blast.run(params.running_dir)

        # TSV generation
        formatter_tsv = BlastFormatter(camel)
        formatter_tsv.update_parameters(output_format='"7 pident sseqid sseq slen qseqid qstart qend"')
        formatter_tsv.add_input_files({'ASN': blast.tool_outputs['ASN']})
        formatter_tsv.run(params.running_dir)

        # Best hit selection
        hit_selector = BestHitSelector(camel)
        hit_selector.add_input_files({'TSV': formatter_tsv.tool_outputs['TSV']})
        hit_selector.add_input_informs({'locus': SnakemakeUtils.load_object(input.INFORMS_locus)})
        hit_selector.run(params.running_dir)

        # Text alignment generation
        formatter_text = BlastFormatter(camel)
        formatter_text.update_parameters(output_format='0', num_alignments=1000)
        formatter_text.add_input_files({'ASN': blast.tool_outputs['ASN']})
        formatter_text.run(params.running_dir)

        # Alignment extraction
        extractor = AlignmentExtractor(camel)
        extractor.add_input_files({'TXT': formatter_text.tool_outputs['TXT'],
                                   'VAL_Hits': hit_selector.tool_outputs['VAL_Hit']})
        extractor.run(params.running_dir)

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
        FASTQ_PE = TRIMMED_READS_PE,
        FASTA = os.path.join(TYPING_WORKING_DIR, '{scheme}', '{locus}', 'fasta-{locus}.io'),
        INFORMS_locus = os.path.join(TYPING_WORKING_DIR, '{scheme}', '{locus}', 'informs-{locus}.io')
    output:
        VAL_Hit = os.path.join(TYPING_WORKING_DIR, '{scheme}', '{locus}', 'hit-srst2.io')
    params:
        running_dir = os.path.join(TYPING_WORKING_DIR, '{scheme}', '{locus}')
    threads:
        4
    run:
        from app.tools.srst2.srst2alleledetector import SRST2AlleleDetector
        detector = SRST2AlleleDetector(camel)
        SnakemakeUtils.add_pickle_inputs(detector, input)
        step = SnakeStep(rule, detector, camel, params.running_dir, config)
        detector.update_parameters(threads=threads, forward_designator='1P', reverse_designator='2P')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(detector, output)

rule combine_hits:
    """
    Combines hits from one scheme into a list. This also ensures that the allele detection is done using blast or
    using SRST2 depending on the config.
    """
    input:
        hits_blast = lambda wildcards: expand(os.path.join(TYPING_WORKING_DIR, '{scheme}', '{locus}', 'hit-blast.io'),
                                              locus=LOCI.get(wildcards.scheme), scheme=wildcards.scheme) if config['detection_method'] == 'fast' else [],
        hits_srst2 = lambda wildcards: expand(os.path.join(TYPING_WORKING_DIR, '{scheme}', '{locus}', 'hit-srst2.io'),
                                              locus=LOCI.get(wildcards.scheme), scheme=wildcards.scheme) if config['detection_method'] == 'normal' else []
    output:
        os.path.join(TYPING_WORKING_DIR, '{scheme}', 'hits-combined.io')
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
        hits=os.path.join(TYPING_WORKING_DIR, '{scheme}', 'hits-combined.io'),
        INFORMS_Loci=lambda wildcards: expand(os.path.join(TYPING_WORKING_DIR, '{scheme}', '{locus}', 'informs-{locus}.io'),
                                              locus=LOCI.get(wildcards.scheme), scheme=wildcards.scheme)
    output:
        TSV=os.path.join(TYPING_WORKING_DIR, '{scheme}', 'tsv-combined.io')
    params:
        running_dir=os.path.join(TYPING_WORKING_DIR, '{scheme}'),
        sample_name=FileSystemHelper.make_valid(config['sample_name']),
        scheme=FileSystemHelper.make_valid('{scheme}')
    run:
        from app.tools.pipelines.sequence_typing.allelecombiner import AlleleCombiner
        list_of_informs = []
        for pickle in input.INFORMS_Loci:
            list_of_informs.append(SnakemakeUtils.load_object(pickle))
        combiner = AlleleCombiner(camel)
        step = SnakeStep(rule, combiner, camel, params.running_dir, config)
        output_path = os.path.join(params.running_dir, 'typing-{}-{}.tsv'.format(params.scheme, params.sample_name))
        combiner.update_parameters(output_filename=output_path)
        combiner.add_input_files({'VAL_Hits': SnakemakeUtils.load_object(input.hits)})
        combiner.add_input_informs({'loci': list_of_informs})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(combiner, output)

rule detect_sequence_type:
    """
    Detects the sequence type based on the detected alleles.
    """
    input:
        VAL_Hits = os.path.join(TYPING_WORKING_DIR, '{scheme}', 'hits-combined.io'),
        TSV = os.path.join(TYPING_WORKING_DIR, '{scheme}', 'tsv-profiles.io')
    output:
        INFORMS = os.path.join(TYPING_WORKING_DIR, '{scheme}', 'informs-ST.io')
    params:
        running_dir = os.path.join(TYPING_WORKING_DIR, '{scheme}')
    run:
        from app.tools.pipelines.sequence_typing.sequencetypedetector import SequenceTypeDetector
        sequence_type_detector = SequenceTypeDetector(camel)
        SnakemakeUtils.add_pickle_inputs(sequence_type_detector, input)
        step = SnakeStep(rule, sequence_type_detector, camel, params.running_dir, config)
        sequence_type_detector.update_parameters(allele_wildcard='N', allele_absent_symbol='0')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sequence_type_detector, output)

rule create_report:
    """
    Creates a report with the sequence typing output.
    """
    input:
        TSV = os.path.join(TYPING_WORKING_DIR, '{scheme}', 'tsv-combined.io'),
        INFORMS_Scheme = os.path.join(TYPING_WORKING_DIR, '{scheme}', 'informs-{scheme}.io'),
        INFORMS_ST = lambda wildcards: os.path.join(TYPING_WORKING_DIR, '{}', 'informs-ST.io').format(
            wildcards.scheme) if has_profiles(wildcards.scheme) else [],
        VAL_Hits = os.path.join(TYPING_WORKING_DIR, '{scheme}', 'hits-combined.io'),
    output:
        VAL_HTML = os.path.join(TYPING_WORKING_DIR, '{scheme}', 'report', 'html.io')
    params:
        running_dir = os.path.join(TYPING_WORKING_DIR, '{scheme}', 'report'),
    run:
        from app.tools.pipelines.sequence_typing.htmlreportertyping import HtmlReporterTyping
        reporter = HtmlReporterTyping(camel)
        if len(input.INFORMS_ST) != 0:
            reporter.add_input_informs({'ST': SnakemakeUtils.load_object(input.INFORMS_ST)})
        SnakemakeUtils.add_pickle_inputs(reporter, input, ['TSV', 'INFORMS_Scheme', 'VAL_Hits'])
        step = SnakeStep(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule create_scheme_summary:
    """
    Create a tabular summary with the sequence typing output on a scheme
    """
    input:
        VAL_Hits=os.path.join(TYPING_WORKING_DIR, '{scheme}', 'hits-combined.io'),
        INFORMS_ST=lambda wildcards: os.path.join(TYPING_WORKING_DIR, wildcards.scheme, 'informs-ST.io') if has_profiles(wildcards.scheme) else []
    output:
        TYPING_SCHEME_SUMMARY
    params:
        scheme_name=lambda wildcards: wildcards.scheme
    run:
        if len(input.INFORMS_ST) == 0:
            st_metadata = []
        else:
            st_metadata = SnakemakeUtils.load_object(input.INFORMS_ST)['sequence_type'].metadata
        hits = SnakemakeUtils.load_object(input.VAL_Hits)
        with open(output[0], 'w') as handle:
            for k, v in st_metadata:
                handle.write(f'{params.scheme_name}-{k}\t{v}')
                handle.write('\n')
            for hit in hits:
                key = '{}-{}'.format(params.scheme_name, hit.value.locus)
                handle.write(f'{key}\t{hit.value.allele_id}')
                handle.write('\n')

rule combine_sequence_typing_reports:
    """
    Combines the reports of all sequence typing schemes.
    """
    input:
        VAL_HTML = expand(os.path.join(TYPING_WORKING_DIR, '{scheme}', 'report', 'html.io'), scheme=__OTHER_ST_DBS)
    output:
        TYPING_REPORT
    params:
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

rule species_confirmation_report:
    input:
        VAL_HTML = expand(os.path.join(TYPING_WORKING_DIR, '{scheme}', 'report', 'html.io'), scheme=__SPECIES_CONFIRM_ST_DBS)
    output:
        SPECIES_CONFIRM_ST_REPORT
    params:
        output_dir = config['output_dir']
    run:
        species_confirm_module = HtmlElement('div')
        species_confirm_module.add_module_header('Species Identification')
        for pickle in input.VAL_HTML:
            report_section = SnakemakeUtils.load_object(pickle)[0].value
            report_section.copy_files(params.output_dir)
            species_confirm_module.add_html_object(report_section)
        SnakemakeUtils.dump_object([ToolIOValue(species_confirm_module)], output[0])

