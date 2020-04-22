from pathlib import Path

from camel.app.camel import Camel
from camel.app.components.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.sequence_typing import OUTPUT_TYPING_HITS

camel = Camel.get_instance()


rule typing_srst2_allele_detection:
    """
    Allele detection using SRST2.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
        INFORMS_scheme = Path(config['working_dir']) / 'typing' / '{scheme}' / 'informs-locus_set.io'
    output:
        VAL_Hit = Path(config['working_dir']) / 'typing' / '{scheme}' / '{locus_type}' / '{locus}' / 'hit-srst2.io',
        INFORMS = Path(config['working_dir']) / 'typing' / '{scheme}' / '{locus_type}' / '{locus}' / 'informs-srst2.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'typing' / wildcards.scheme / wildcards.locus_type / wildcards.locus,
        locus_name = lambda wildcards: wildcards.locus,
        scheme_dir = lambda wildcards: Path(SCHEME_DATA[wildcards.scheme]['path']),
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE',
        srst2_options = config.get('srst2')
    threads: 4
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.srst2.srst2alleledetector import SRST2AlleleDetector

        # Get metadata
        scheme_informs = SnakemakeUtils.load_object(input.INFORMS_scheme)
        locus_informs = scheme_informs['loci'].metadata_by_locus_name[params.locus_name]

        detector = SRST2AlleleDetector(camel)
        fastq_input = SnakePipelineUtils.extracts_fq_input(
            input.IO, key_pe='FASTQ_PE', key_se='FASTQ_SE', read_type=params.read_type)
        detector.add_input_files(fastq_input)
        detector.add_input_files({'FASTA': [ToolIOFile(str(params.scheme_dir / locus_informs['fasta_path']))]})
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

rule typing_srst2_combine_hits:
    """
    Combines the separate SRST2 hits for a scheme into a single IO object.
    """
    input:
        IO_hits = lambda wildcards: [str(rules.typing_srst2_allele_detection.output.VAL_Hit).format(locus=locus, scheme=wildcards.scheme, locus_type='DNA') for locus in loci_by_scheme_by_type[wildcards.scheme]['DNA']],
    output:
        IO = Path(config['working_dir']) / str(OUTPUT_TYPING_HITS).format(locus_type='DNA', scheme='{scheme}', detection_method='srst2')
    run:
        all_hits = []
        for io in input.IO_hits:
            all_hits.extend(SnakemakeUtils.load_object(io))
        SnakemakeUtils.dump_object(all_hits, str(output.IO))
