import logging

import os

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.assembly_spades import OUTPUT_ASSEMBLY_FASTA
from camel.resources.snakefile.sequence_typing import OUTPUT_TYPING_HITS

camel = Camel.get_instance()


rule Typing_blast_allele_detection:
    """
    Allele detection using blastn (DNA) or blastx (peptide).
    """
    input:
        FASTA=os.path.join(config['working_dir'], OUTPUT_ASSEMBLY_FASTA),
        INFORMS_scheme=os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-locus_set.io')
    output:
        VAL_Hit=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', '{locus}', 'hit-blast.io')
    params:
        working_dir=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', '{locus}'),
        locus_type=lambda wildcards: wildcards.locus_type,
        locus_name=lambda wildcards: wildcards.locus,
        scheme_dir=lambda wildcards: SCHEMES[wildcards.scheme]
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
        elif params.locus_type == 'peptide':
            blast = Blastx(camel)
        else:
            raise ValueError(f"Invalid locus type: {wildcards.locus_type}")
        SnakemakeUtils.add_pickle_input(blast, 'FASTA', input.FASTA)
        blast.add_input_files({'DB_BLAST': [ToolIOFile(os.path.join(params.scheme_dir, locus_informs['fasta_path']))]})
        blast.update_parameters(threads=threads)
        blast.run(params.working_dir)

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

rule Typing_blast_combine_hits:
    """
    Combines the hits for the blast based detection.
    """
    input:
        input_nucl=lambda wildcards: expand(os.path.join(config['working_dir'], 'typing', wildcards.scheme, 'DNA', '{locus}', 'hit-blast.io'), locus=loci_by_scheme_by_type[wildcards.scheme]['DNA']),
        cleanup_nucl=lambda wildcards: expand(os.path.join(config['working_dir'], 'typing', wildcards.scheme, 'DNA', '{locus}', '.cleanup'), locus=loci_by_scheme_by_type[wildcards.scheme]['DNA']),
        input_pept=lambda wildcards: expand(os.path.join(config['working_dir'], 'typing', wildcards.scheme, 'peptide', '{locus}', 'hit-blast.io'), locus=loci_by_scheme_by_type[wildcards.scheme]['peptide']),
        cleanup_pept=lambda wildcards: expand(os.path.join(config['working_dir'], 'typing', wildcards.scheme, 'peptide', '{locus}', '.cleanup'), locus=loci_by_scheme_by_type[wildcards.scheme]['peptide'])
    output:
        hits_nucl=os.path.join(config['working_dir'], OUTPUT_TYPING_HITS.format(scheme='{scheme}', detection_method='blast', locus_type='DNA')),
        hits_pept=os.path.join(config['working_dir'], OUTPUT_TYPING_HITS.format(scheme='{scheme}', detection_method='blast', locus_type='peptide'))
    run:
        for key in 'nucl', 'pept':
            list_of_hits = []
            for pickle in input.get(f'input_{key}'):
                list_of_hits.append(SnakemakeUtils.load_object(pickle)[0])
            SnakemakeUtils.dump_object(list_of_hits, output.get(f'hits_{key}'))

rule Typing_blastn_cleanup:
    """
    Cleans up the directories created by the blast allele calling.
    """
    input:
        VAL_Hit=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', '{locus}', 'hit-blast.io')
    output:
        flag=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', '{locus}', '.cleanup')
    priority: 5
    run:
        import humanize
        dir_analysis = os.path.dirname(input.VAL_Hit)
        bytes_cleared = 0
        removed_extensions = ['.asn', '.bam', '.pileup']
        for f in os.listdir(dir_analysis):
            if os.path.splitext(f)[-1] not in removed_extensions:
                continue
            full_path = os.path.join(dir_analysis, f)
            bytes_cleared += os.path.getsize(full_path)
            os.remove(full_path)
        logging.info('{} cleaned'.format(humanize.naturalsize(bytes_cleared)))
        with open(output.flag, 'w') as handle_out:
            handle_out.write('done')
