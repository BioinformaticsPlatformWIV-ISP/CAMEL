import os

from camel.app.components.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_READS_PE
from camel.resources.snakefile.read_trimming_iontorrent import OUTPUT_TRIMMING_IT_READS
from camel.resources.snakefile.sequence_typing import OUTPUT_TYPING_HITS


rule Typing_srst2_select_input:
    """
    Selects the input for the SRST2 sequence typing. 
    """
    input:
        FASTQ_PE=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_PE) if config.get('read_type', 'illumina') == 'illumina' else [],
        FASTQ_SE=os.path.join(config['working_dir'], OUTPUT_TRIMMING_IT_READS) if config.get('read_type', 'illumina') == 'iontorrent' else []
    output:
        FASTQ=os.path.join(config['working_dir'], 'typing', 'input-fastq.io')
    params:
        read_type=config.get('read_type', 'illumina')
    run:
        input_dict = {}
        for key in input.keys():
            if len(input[key]) > 0:
                input_dict[key] = SnakemakeUtils.load_object(input[key])
        SnakemakeUtils.dump_object(input_dict, output[0])

rule Typing_srst2_allele_detection:
    """
    Allele detection using SRST2.
    """
    input:
        FASTQ = os.path.join(config['working_dir'], 'typing', 'input-fastq.io'),
        INFORMS_scheme = os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-locus_set.io')
    output:
        VAL_Hit = os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', '{locus}', 'hit-srst2.io')
    params:
        running_dir = os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', '{locus}'),
        locus_name = lambda wildcards: wildcards.locus,
        scheme_dir = lambda wildcards: SCHEMES[wildcards.scheme],
        srst2_options = config.get('srst2')
    threads: 4
    run:
        from camel.app.tools.srst2.srst2alleledetector import SRST2AlleleDetector

        # Get metadata
        scheme_informs = SnakemakeUtils.load_object(input.INFORMS_scheme)
        locus_informs = scheme_informs['loci'].metadata_by_locus_name[params.locus_name]

        detector = SRST2AlleleDetector(camel)
        fastq_input = SnakemakeUtils.load_object(input.FASTQ)
        detector.add_input_files(fastq_input)
        detector.add_input_files({'FASTA': [ToolIOFile(os.path.join(params.scheme_dir, locus_informs['fasta_path']))]})
        detector.add_input_informs({'locus': locus_informs})
        if (params.srst2_options is not None) and ('max_unaligned_overlap' in params.srst2_options):
            detector.update_parameters(max_unaligned_overlap=params.srst2_options['max_unaligned_overlap'])
        step = Step(rule, detector, camel, params.running_dir, config)
        if 'FASTQ_PE' in fastq_input:
            fwd_read_path = fastq_input['FASTQ_PE'][0].path
            fwd_designator, rev_designator = SequenceTypingUtils.determine_read_status(fwd_read_path)
            detector.update_parameters(forward_designator=fwd_designator, reverse_designator=rev_designator)
        detector.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(detector, output)

rule Typing_srst2_combine_hits:
    """
    Combines the separate SRST2 hits into a single IO object.
    """
    input:
        input_nucl=lambda wildcards: expand(os.path.join(config['working_dir'], 'typing', wildcards.scheme, 'DNA', '{locus}', 'hit-srst2.io'), locus=loci_by_scheme_by_type[wildcards.scheme]['DNA'])
    output:
        os.path.join(config['working_dir'], OUTPUT_TYPING_HITS.format(locus_type='DNA', scheme='{scheme}', detection_method='srst2'))
    run:
        list_of_hits = []
        for pickle in input:
            list_of_hits += SnakemakeUtils.load_object(pickle)
        SnakemakeUtils.dump_object(list_of_hits, output[0])
