from pathlib import Path

from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.sequence_typing import OUTPUT_TYPING_HITS



rule typing_kma_allele_detection:
    """
    Allele detection using KMA.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
        INFORMS_scheme = Path(config['working_dir']) / 'typing' / '{scheme}' / 'informs-locus_set.io'
    output:
        VAL_hit = Path(config['working_dir']) / 'typing' / '{scheme}' / '{locus_type}' / '{locus}' / 'hit-kma.io',
        INFORMS = Path(config['working_dir']) / 'typing' / '{scheme}' / '{locus_type}' / '{locus}' / 'informs-kma.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'typing' / wildcards.scheme / wildcards.locus_type / wildcards.locus,
        locus_name = lambda wildcards: wildcards.locus,
        scheme_dir = lambda wildcards: SCHEME_DATA[wildcards.scheme]['path']
    threads: 4
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.kma.kma import KMA
        from camel.app.tools.kma.kmatypinghitextractor import KMATypingHitExtractor

        # Get metadata
        scheme_informs = SnakemakeUtils.load_object(input.INFORMS_scheme)
        locus_informs = scheme_informs['loci'].metadata_by_locus_name[params.locus_name]
        dir_kma = (Path(str(params.scheme_dir)) / locus_informs['fasta_path']).parent / 'kma'
        try:
            db_path = next(dir_kma.glob('*.name'))
        except StopIteration:
            raise FileNotFoundError(f"KMA database for locus '{params.locus_name}' ({params.scheme_dir}) not found")

        # Launch KMA
        kma = KMA(camel)
        fastq_input = SnakePipelineUtils.extracts_fq_input(input.IO, key_pe='FASTQ_PE')
        kma.add_input_files(fastq_input)
        kma.add_input_files({'DB': [ToolIOValue(str(db_path.parent / db_path.stem))]})
        step = Step(rule, kma, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_object(kma.informs, output.INFORMS)

        # Extract the best hit
        kma_extractor = KMATypingHitExtractor(camel)
        kma_extractor.add_input_files({'TSV': kma.tool_outputs['TSV']})
        kma_extractor.add_input_informs({'locus': locus_informs})
        kma_extractor.run(str(params.running_dir))
        SnakemakeUtils.dump_tool_output(kma_extractor, 'VAL_hit', output.VAL_hit)

rule typing_kma_combine_hits:
    """
    Combines the separate SRST2 hits into a single IO object.
    """
    input:
        VAL_HIT_NUCL = lambda wildcards: expand(str(rules.typing_kma_allele_detection.output.VAL_hit), locus=loci_by_scheme_by_type[wildcards.scheme]['DNA'], scheme='{scheme}', locus_type='DNA')
    output:
        VAL_HITS = Path(config['working_dir']) / str(OUTPUT_TYPING_HITS).format(locus_type='DNA', scheme='{scheme}', detection_method='kma')
    run:
        list_of_hits = []
        for pickle in input.VAL_HIT_NUCL:
            list_of_hits += SnakemakeUtils.load_object(pickle)
        SnakemakeUtils.dump_object(list_of_hits, output.VAL_HITS)
