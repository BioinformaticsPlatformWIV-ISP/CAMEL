import os

from camel.app.components.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_READS_PE
from camel.resources.snakefile.sequence_typing import OUTPUT_TYPING_HITS


rule Typing_srst2_select_input:
    """
    Selects the input file for SRST2.
    """
    input:
        FASTQ_PE=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_PE) if config.get('read_type', 'illumina') == 'illumina' else []
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
        FASTQ=os.path.join(config['working_dir'], 'typing', 'input-fastq.io'),
        FASTA=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', '{locus}', 'fasta.io'),
        INFORMS_locus=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', '{locus}', 'informs.io')
    output:
        VAL_Hit=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', '{locus}', 'hit-srst2.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', '{locus}')
    threads: 4
    run:
        from camel.app.tools.srst2.srst2alleledetector import SRST2AlleleDetector
        detector = SRST2AlleleDetector(camel)
        SnakemakeUtils.add_pickle_inputs(detector, input, keys=['FASTA', 'INFORMS_locus'])
        fastq_input = SnakemakeUtils.load_object(input.FASTQ)
        detector.add_input_files(fastq_input)
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
