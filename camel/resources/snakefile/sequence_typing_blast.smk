from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.assembly_spades import OUTPUT_ASSEMBLY_FASTA
from camel.resources.snakefile.sequence_typing import OUTPUT_TYPING_HITS

camel = Camel.get_instance()


rule typing_blast_allele_detection:
    """
    Allele detection using blastn (DNA) or blastx (peptide).
    """
    input:
        FASTA = Path(config['working_dir']) / OUTPUT_ASSEMBLY_FASTA,
        INFORMS_scheme = Path(config['working_dir']) / 'typing' / '{scheme}' / 'informs-locus_set.io'
    output:
        VAL_Hit = Path(config['working_dir']) / 'typing' / '{scheme}' / '{locus_type}' / '{locus}' / 'hit-blast.io',
        INFORMS = Path(config['working_dir']) / 'typing' / '{scheme}' / '{locus_type}' / '{locus}' / 'informs-blast.io'
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / 'typing' / wildcards.scheme / wildcards.locus_type / wildcards.locus,
        locus_type = lambda wildcards: wildcards.locus_type,
        locus_name = lambda wildcards: wildcards.locus,
        scheme_dir = lambda wildcards: Path(SCHEME_DATA[wildcards.scheme]['path']),
        blastn_task = lambda wildcards: SCHEME_DATA[wildcards.scheme].get('blastn_task', 'megablast')
    threads: 1
    run:
        from camel.app.tools.blast.blastformatter import BlastFormatter
        from camel.app.tools.blast.blastn import Blastn
        from camel.app.tools.blast.blastx import Blastx
        from camel.app.tools.pipelines.sequence_typing.besthitselector import BestHitSelector
        from camel.app.tools.pipelines.sequence_typing.alignmentextractor import AlignmentExtractor
        from camel.app.components.blast.blasthitstatistics import BLASTN_OUTPUT_FORMAT

        # Get metadata
        scheme_informs = SnakemakeUtils.load_object(input.INFORMS_scheme)
        locus_informs = scheme_informs['loci'].metadata_by_locus_name[params.locus_name]

        # Blast alignment
        if params.locus_type == 'DNA':
            blast = Blastn(camel)
            blast.update_parameters(task=params.blastn_task)
        elif params.locus_type == 'peptide':
            blast = Blastx(camel)
        else:
            raise ValueError(f"Invalid locus type: {wildcards.locus_type}")
        SnakemakeUtils.add_pickle_input(blast, 'FASTA', input.FASTA)
        blast.add_input_files({'DB_BLAST': [ToolIOFile(str(params.scheme_dir / locus_informs['fasta_path']))]})
        blast.update_parameters(threads=threads)
        blast.run(str(params.working_dir))
        SnakemakeUtils.dump_object(blast.informs, output.INFORMS)

        # TSV generation
        formatter_tsv = BlastFormatter(camel)
        formatter_tsv.update_parameters(output_format=BLASTN_OUTPUT_FORMAT)
        formatter_tsv.add_input_files({'ASN': blast.tool_outputs['ASN']})
        formatter_tsv.run(params.working_dir)

        # Best hit selection
        hit_selector = BestHitSelector(camel)
        hit_selector.add_input_files({'TSV': formatter_tsv.tool_outputs['TSV']})
        hit_selector.add_input_informs({'locus': locus_informs})
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

rule typing_blast_combine_hits:
    """
    Combines the hits for the blast based detection.
    """
    input:
        input_nucl=lambda wildcards: expand(str(Path(config['working_dir']) / 'typing' / wildcards.scheme / 'DNA' / '{locus}' / 'hit-blast.io'), locus=loci_by_scheme_by_type[wildcards.scheme]['DNA']),
        input_pept=lambda wildcards: expand(str(Path(config['working_dir']) / 'typing' / wildcards.scheme / 'peptide' / '{locus}' / 'hit-blast.io'), locus=loci_by_scheme_by_type[wildcards.scheme]['peptide']),
    output:
        hits_nucl = Path(config['working_dir']) / str(OUTPUT_TYPING_HITS).format(scheme='{scheme}', detection_method='blast', locus_type='DNA'),
        hits_pept = Path(config['working_dir']) / str(OUTPUT_TYPING_HITS).format(scheme='{scheme}', detection_method='blast', locus_type='peptide')
    run:
        for key in 'nucl', 'pept':
            list_of_hits = []
            for pickle in input.get(f'input_{key}'):
                list_of_hits.append(SnakemakeUtils.load_object(pickle)[0])
            SnakemakeUtils.dump_object(list_of_hits, output.get(f'hits_{key}'))
