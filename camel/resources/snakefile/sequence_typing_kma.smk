import os
from pathlib import Path

from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.sequence_typing import OUTPUT_TYPING_HITS



rule Typing_KMA_allele_detection:
    """
    Allele detection using KMA.
    """
    input:
        FASTQ_PE = os.path.join(config['working_dir'], 'typing', 'input-fastq.io'),
        INFORMS_scheme = os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-locus_set.io')
    output:
        VAL_hit = os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', '{locus}', 'hit-kma.io')
    params:
        running_dir = os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', '{locus}'),
        locus_name = lambda wildcards: wildcards.locus,
        scheme_dir = lambda wildcards: SCHEMES[wildcards.scheme]
    threads: 4
    run:
        from camel.app.tools.kma.kma import KMA
        from camel.app.tools.kma.kmatypinghitextractor import KMATypingHitExtractor

        # Get metadata
        scheme_informs = SnakemakeUtils.load_object(input.INFORMS_scheme)
        locus_informs = scheme_informs['loci'].metadata_by_locus_name[params.locus_name]
        locus_name_valid = locus_informs['name_valid']
        db_path = Path(params.scheme_dir) / locus_informs['kma_path']
        if not db_path.exists():
            raise FileNotFoundError(f"KMA database for '{locus_name}' not found")

        # Launch KMA
        kma = KMA(camel)
        kma.add_input_files(SnakemakeUtils.load_object(input.FASTQ_PE))
        kma.add_input_files({'DB': [ToolIOValue(str(db_path))]})
        step = Step(rule, kma, camel, params.running_dir, config)
        step.run_step()

        # Extract the best hit
        kma_extractor = KMATypingHitExtractor(camel)
        kma_extractor.add_input_files({'TSV': kma.tool_outputs['TSV']})
        kma_extractor.add_input_informs({'locus': locus_informs})
        kma_extractor.run(params.running_dir)
        SnakemakeUtils.dump_tool_output(kma_extractor, 'VAL_hit', output.VAL_hit)

rule Typing_KMA_combine_hits:
    """
    Combines the separate SRST2 hits into a single IO object.
    """
    input:
        input_nucl=lambda wildcards: expand(os.path.join(config['working_dir'], 'typing', wildcards.scheme, 'DNA', '{locus}', 'hit-kma.io'), locus=loci_by_scheme_by_type[wildcards.scheme]['DNA'])
    output:
        os.path.join(config['working_dir'], OUTPUT_TYPING_HITS.format(locus_type='DNA', scheme='{scheme}', detection_method='kma'))
    run:
        list_of_hits = []
        for pickle in input:
            list_of_hits += SnakemakeUtils.load_object(pickle)
        SnakemakeUtils.dump_object(list_of_hits, output[0])
