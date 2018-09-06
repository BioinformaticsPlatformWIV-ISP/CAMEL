"""
This Snakefile can be used to perform sequence typing.
There are two approaches:
- Fast: blastn for nucleotide loci and blastx for protein loci.
- Normal: SRST2 for nucleotide loci and blastx for protein loci.

The sequence typing databases need to be added to the config file under the 'sequence_typing' section.
Each database need to be specified as 'Key': 'Path to DB folder'.
"""
from typing import List

import os
import json

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
# TODO: Re-use these imports when those Snakefiles are available
# from camel.resources.snakefile.assembly_spades import OUTPUT_ASSEMBLY_FASTA
# from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_READS_PE
OUTPUT_ASSEMBLY_FASTA = os.path.join(config['working_dir'], 'sequence_typing', 'input', 'fasta.io')
OUTPUT_READ_TRIMMING_READS_PE = os.path.join(config['working_dir'], 'sequence_typing', 'input', 'fastq.io')
# ------
from camel.resources.snakefile.sequence_typing import OUTPUT_TYPING_REPORT, OUTPUT_TYPING_REPORT_EMPTY, OUTPUT_TYPING_SUMMARY




camel = Camel()


def get_loci(folder: str) -> List[str]:
    """
    Returns the loci from the given database.
    :param folder: Database folder.
    :return: Loci from database
    """
    return sorted([d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d)) and not d.startswith('.')])


def has_profiles(scheme_key: str) -> bool:
    """
    Returns True if the scheme has profiles.
    :param scheme_key: Scheme
    :return: True if the scheme has profiles
    """
    if scheme_key not in SCHEMES:
        raise ValueError(f"Scheme '{scheme_key}' not found")
    folder = SCHEMES.get(scheme_key)
    return os.path.isfile(os.path.join(folder, 'profiles.tsv'))


def get_locus_type(scheme: str, locus: str) -> str:
    """
    Returns the type of locus ('DNA', 'peptide').
    :param scheme: Scheme (e.g. 'MLST')
    :param locus: Locus (e.g. 'abcZ')
    :return: Locus type
    """
    locus_directory = os.path.join(SCHEMES[scheme], locus)
    locus_metadata_file = os.path.join(locus_directory, 'locus_metadata.txt')
    if not os.path.isfile(locus_metadata_file):
        raise FileNotFoundError("No metadata found in '{}'".format(locus_directory))
    with open(locus_metadata_file) as handle:
        try:
            return json.load(handle)['type']
        except KeyError:
            raise ValueError("Metadata file does not contain locus type ({})".format(locus_directory))


LOCI = {key: get_loci(path) for key, path in config['sequence_typing'].items()}
SCHEMES = config['sequence_typing']

# Supported detection methods:
# - blast: BLASTN / BLASTX allele detection
# - srst2: SRST2 allele detection
DETECTION_METHOD = config['detection_method']


rule Sequence_typing_extract_locus_set_info:
    """
    Extracts the scheme metadata.
    """
    input:
        lambda wildcards: SCHEMES[wildcards.scheme]
    output:
        INFORMS=os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-locus_set.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'typing', '{scheme}')
    run:
        from camel.app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
        locus_set_manager = LocusSetManager(camel)
        locus_set_manager.add_input_files({'DIR': [ToolIODirectory(input[0])]})
        step = Step(rule, locus_set_manager, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_object(locus_set_manager.informs, output.INFORMS)

rule Sequence_typing_extract_locus_info:
    """
    Extracts the metadata for a single locus.
    """
    input:
        lambda wildcards: os.path.join(SCHEMES.get(wildcards.scheme), wildcards.locus)
    output:
        FASTA=os.path.join(config['working_dir'], 'typing', '{scheme}', 'fasta-{locus}.io'),
        INFORMS=os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-{locus}.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'typing', '{scheme}')
    run:
        from camel.app.tools.pipelines.sequence_typing.locusmanager import LocusManager
        locus_manager = LocusManager(camel)
        locus_manager.add_input_files({'DIR': [ToolIODirectory(input[0])]})
        step = Step(rule, locus_manager, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(locus_manager, output)

rule Sequence_typing_pickle_dump_sequence_type_definitions:
    """
    Retrieves the sequence type definitions and converts them to CAMEL IO pickles.
    """
    input:
        lambda wildcards: SCHEMES[wildcards.scheme]
    output:
        os.path.join(config['working_dir'], 'typing', '{scheme}', 'tsv-profiles.io')
    run:
        SnakemakeUtils.dump_object([ToolIOFile(os.path.join(input[0], 'profiles.tsv'))], output[0])

rule Sequence_typing_blast_allele_detection:
    """
    Allele detection using blastn.
    """
    input:
        FASTA=os.path.join(config['working_dir'], OUTPUT_ASSEMBLY_FASTA),
        DB_BLAST=os.path.join(config['working_dir'], 'typing', '{scheme}', 'fasta-{locus}.io'),
        INFORMS_locus=os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-{locus}.io')
    output:
        VAL_Hit=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus}', 'hit-blast.io')
    params:
        working_dir=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus}')
    threads: 1
    run:
        from camel.app.tools.blast.blastformatter import BlastFormatter
        from camel.app.tools.blast.blastn import Blastn
        from camel.app.tools.blast.blastx import Blastx
        from camel.app.tools.pipelines.sequence_typing.besthitselector import BestHitSelector
        from camel.app.tools.pipelines.sequence_typing.alignmentextractor import AlignmentExtractor

        # Blast alignment
        locus_type = get_locus_type(wildcards.scheme, wildcards.locus)
        if locus_type == 'DNA':
            blast = Blastn(camel)
        elif locus_type == 'peptide':
            blast = Blastx(camel)
        else:
            raise ValueError(f"Invalid locus type: {locus_type}")
        SnakemakeUtils.add_pickle_input(blast, 'DB_BLAST', input.DB_BLAST)
        SnakemakeUtils.add_pickle_input(blast, 'FASTA', input.FASTA)
        blast.update_parameters(threads=threads)
        blast.run(params.working_dir)

        # TSV generation
        formatter_tsv = BlastFormatter(camel)
        formatter_tsv.update_parameters(output_format='"7 pident sseqid sseq slen qseqid qstart qend"')
        formatter_tsv.add_input_files({'ASN': blast.tool_outputs['ASN']})
        formatter_tsv.run(params.working_dir)

        # Best hit selection
        hit_selector = BestHitSelector(camel)
        hit_selector.add_input_files({'TSV': formatter_tsv.tool_outputs['TSV']})
        hit_selector.add_input_informs({'locus': SnakemakeUtils.load_object(input.INFORMS_locus)})
        hit_selector.run(params.working_dir)

        # Text alignment generation
        formatter_text = BlastFormatter(camel)
        formatter_text.update_parameters(output_format='0', num_alignments=1000)
        formatter_text.add_input_files({'ASN': blast.tool_outputs['ASN']})
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

rule Sequence_typing_combine_blast_hits:
    """
    Combines the separate BLAST hits into a single IO object.
    """
    input:
        lambda wildcards: expand(os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus}', 'hit-blast.io'),
                                 locus=LOCI.get(wildcards.scheme),
                                 scheme=wildcards.scheme)
    output:
        os.path.join(config['working_dir'], 'typing', '{scheme}', 'blast', 'all-hits.io')
    run:
        list_of_hits = []
        for pickle in input:
            list_of_hits += SnakemakeUtils.load_object(pickle)
        SnakemakeUtils.dump_object(list_of_hits, output[0])

rule Sequence_typing_srst2_allele_detection:
    """
    Allele detection using SRST2.
    """
    input:
        FASTQ_PE=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_PE),
        FASTA=os.path.join(config['working_dir'], 'typing', '{scheme}', 'fasta-{locus}.io'),
        INFORMS_locus=os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-{locus}.io')
    output:
        VAL_Hit=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus}', 'hit-srst2.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus}')
    threads: 4
    run:
        from camel.app.tools.srst2.srst2alleledetector import SRST2AlleleDetector
        detector = SRST2AlleleDetector(camel)
        SnakemakeUtils.add_pickle_inputs(detector, input)
        step = Step(rule, detector, camel, params.running_dir, config)
        detector.update_parameters(threads=threads, forward_designator='1P', reverse_designator='2P')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(detector, output)

rule Sequence_typing_combine_srst2_hits:
    """
    Combines the separate SRST2 hits into a single IO object.
    """
    input:
        lambda wildcards: expand(os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus}', 'hit-srst2.io'),
                                 locus=LOCI.get(wildcards.scheme),
                                 scheme=wildcards.scheme)
    output:
        os.path.join(config['working_dir'], 'typing', '{scheme}', 'srst2', 'all-hits.io')
    run:
        list_of_hits = []
        for pickle in input:
            list_of_hits += SnakemakeUtils.load_object(pickle)
        SnakemakeUtils.dump_object(list_of_hits, output[0])

rule Sequence_typing_KMA_allele_detection:
    """
    Allele detection using KMA.
    """
    input:
        FASTQ=OUTPUT_READ_TRIMMING_READS_PE,
        INFORMS_locus=os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-{locus}.io')
    output:
        VAL_hit=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus}', 'hit-kma.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus}')
    run:
        from camel.app.tools.kma.kma import KMA
        kma = KMA(camel)
        SnakemakeUtils.add_pickle_input(kma, 'FASTQ', input.FASTQ)
        kma_path = os.path.join('/scratch/bebog/kma/', wildcards.scheme,
                                wildcards.locus, wildcards.locus)
        kma.add_input_files({'DB_KMA': [ToolIOValue(kma_path)]})
        step = Step(rule, kma, camel, params.running_dir, config)
        step.run_step()

        from camel.app.tools.kma.kmahitextractor import KMAHitExtractor
        extractor = KMAHitExtractor(camel)
        extractor.add_input_informs({'locus': SnakemakeUtils.load_object(input.INFORMS_locus)})
        extractor.add_input_files({'TSV': kma.tool_outputs['TSV']})
        step = Step(rule, extractor, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(extractor, output)

rule Sequence_typing_Combine_KMA_hits:
    """
    Combines the separate SRST2 hits into a single IO object.
    """
    input:
        lambda wildcards: expand(os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus}', 'hit-kma.io'),
                                 locus=LOCI.get(wildcards.scheme),
                                 scheme=wildcards.scheme)
    output:
        os.path.join(config['working_dir'], 'typing', '{scheme}', 'kma', 'all-hits.io')
    run:
        list_of_hits = []
        for pickle in input:
            list_of_hits += SnakemakeUtils.load_object(pickle)
        SnakemakeUtils.dump_object(list_of_hits, output[0])

rule Sequence_typing_combine_loci:
    """
    Combines the output of all loci of a scheme into a tabular output file.
    """
    input:
        hits=os.path.join(config['working_dir'], 'typing', '{scheme}', DETECTION_METHOD, 'all-hits.io'),
        INFORMS_Loci=lambda wildcards: expand(os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-{locus}.io'),
                                              locus=LOCI.get(wildcards.scheme), scheme=wildcards.scheme)
    output:
        TSV=os.path.join(config['working_dir'], 'typing', '{scheme}', 'tsv-combined.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'typing', '{scheme}'),
        sample_name=FileSystemHelper.make_valid(config['sample_name']),
        scheme=lambda wildcards: FileSystemHelper.make_valid(wildcards.scheme)
    run:
        from camel.app.tools.pipelines.sequence_typing.allelecombiner import AlleleCombiner
        list_of_informs = []
        for pickle in input.INFORMS_Loci:
            list_of_informs.append(SnakemakeUtils.load_object(pickle))
        combiner = AlleleCombiner(camel)
        step = Step(rule, combiner, camel, params.running_dir, config)
        output_path = os.path.join(params.running_dir, 'typing-{}-{}.tsv'.format(params.scheme, params.sample_name))
        combiner.update_parameters(output_filename=output_path)
        combiner.add_input_files({'VAL_Hits': SnakemakeUtils.load_object(input.hits)})
        combiner.add_input_informs({'loci': list_of_informs})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(combiner, output)

rule Sequence_typing_detect_sequence_type:
    """
    Detects the sequence type based on the detected alleles.
    """
    input:
        VAL_Hits=os.path.join(config['working_dir'], 'typing', '{scheme}', DETECTION_METHOD, 'all-hits.io'),
        TSV=os.path.join(config['working_dir'], 'typing', '{scheme}', 'tsv-profiles.io')
    output:
        INFORMS=os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-st.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'typing', '{scheme}')
    run:
        from camel.app.tools.pipelines.sequence_typing.sequencetypedetector import SequenceTypeDetector
        sequence_type_detector = SequenceTypeDetector(camel)
        SnakemakeUtils.add_pickle_inputs(sequence_type_detector, input)
        step = Step(rule, sequence_type_detector, camel, params.running_dir, config)
        sequence_type_detector.update_parameters(allele_wildcard='N', allele_absent_symbol='0')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sequence_type_detector, output)

rule Sequence_typing_create_report:
    """
    Creates a report with the sequence typing output.
    """
    input:
        TSV=os.path.join(config['working_dir'], 'typing', '{scheme}', 'tsv-combined.io'),
        INFORMS_scheme=os.path.join(config['working_dir'], 'typing' , '{scheme}', 'informs-locus_set.io'),
        INFORMS_ST=lambda wildcards: os.path.join(config['working_dir'], 'typing',  wildcards.scheme, 'informs-st.io') if has_profiles(wildcards.scheme) else [],
        VAL_Hits=os.path.join(config['working_dir'], 'typing', '{scheme}', DETECTION_METHOD, 'all-hits.io')
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_TYPING_REPORT)
    params:
        running_dir=os.path.join(config['working_dir'], 'typing', '{scheme}'),
        sample_name=config['sample_name']
    run:
        from camel.app.tools.pipelines.sequence_typing.htmlreportertyping import HtmlReporterTyping
        reporter = HtmlReporterTyping(camel)
        if len(input.INFORMS_ST) != 0:
            reporter.add_input_informs({'ST': SnakemakeUtils.load_object(input.INFORMS_ST)})
        SnakemakeUtils.add_pickle_inputs(reporter, input, ['TSV', 'INFORMS_scheme', 'VAL_Hits'])
        reporter.add_input_files({'VAL_SAMPLE': [ToolIOValue(params.sample_name)]})
        step = Step(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule Sequence_typing_create_report_empty:
    """
    Creates an empty sequence typing report when the analysis is disabled.
    """
    input:
        INFORMS_Scheme=os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-locus_set.io')
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_TYPING_REPORT_EMPTY)
    run:
        informs = SnakemakeUtils.load_object(input.INFORMS_Scheme)
        section = HtmlReportSection(informs['title'], 3)
        section.add_paragraph('Analysis disabled')
        SnakemakeUtils.dump_object([ToolIOValue(section)], output[0])

rule Sequence_typing_dump_summary_info:
    """
    Dumps the summary information in tabular format.
    """
    input:
        VAL_Hits=os.path.join(config['working_dir'], 'typing', '{scheme}', DETECTION_METHOD, 'all-hits.io'),
        INFORMS_ST=lambda wildcards: os.path.join(config['working_dir'], 'typing', wildcards.scheme, 'informs-st.io') if has_profiles(wildcards.scheme) else [],
    output:
        os.path.join(config['working_dir'], OUTPUT_TYPING_SUMMARY)
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
