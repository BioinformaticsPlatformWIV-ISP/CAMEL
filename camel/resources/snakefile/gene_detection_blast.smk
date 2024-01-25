from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import gene_detection

rule gene_detection_blast_blastn:
    """
    Aligns sequences agains the database with blastn.
    """
    input:
        FASTA = Path(config['working_dir']) / gene_detection.INPUT_GENE_DETECTION_FASTA,
        DB_BLAST = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'db_manager' / 'fasta-clust.io'
    output:
        ASN = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'blastn' / 'asn.io',
        INFORMS = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'blastn' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'blastn',
        task = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('blastn', {}).get('task', 'megablast'),
        blast_reads= lambda wildcards: config['gene_detection'][wildcards.db].get('params',{}).get('blastn',{}).get('blast_reads', False),
    run:
        from camel.app.tools.blast.blastn import Blastn
        blastn = Blastn(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(blastn, input)
        step = Step(str(rule), blastn, Camel.get_instance(), Path(str(params.running_dir)), wildcards)
        if params.blast_reads:
            blastn.update_parameters(threads=1,task=str(params.task),max_target_seqs=1) #report only 1 hit per read
        else:
            blastn.update_parameters(threads=1, task=str(params.task), max_target_seqs=20000)
        step.run_step()
        blastn.informs['Task'] = params.task
        SnakemakeUtils.dump_tool_outputs(blastn, output)

rule gene_detection_blast_tsv_generation:
    """
    Generates tabular output format to extract hit statistics.
    """
    input:
        ASN = rules.gene_detection_blast_blastn.output.ASN
    output:
        TSV = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'tsv_generation' / 'tsv.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'tsv_generation'
    run:
        from camel.app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(blast_formatter, input)
        step = Step(str(rule), blast_formatter, Camel.get_instance(), Path(str(params.running_dir)), wildcards)
        blast_formatter.update_parameters(output_format='"7 pident sseqid sseq slen qseqid qstart qend score"')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blast_formatter, output)

rule gene_detection_blast_hit_filtering:
    """
    Filters hits based on percent identity and query coverage.
    Extracts the hit information based on the database metadata.
    """
    input:
        TSV = rules.gene_detection_blast_tsv_generation.output.TSV,
        INFORMS_blastn = rules.gene_detection_blast_blastn.output.INFORMS
    output:
        VAL_Hits = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'hit_filtering' / 'blast-hits.io',
        INFORMS = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS_METHOD).format(db='{db}', method='blast')
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'hit_filtering',
        min_percent_identity = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('blastn', {}).get('min_percent_identity', 90),
        min_coverage = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('blastn', {}).get('min_coverage', 60),
        filtering_method = lambda wildcards: config['gene_detection'][wildcards.db].get('params',{}).get('blastn',{}).get('filtering_method', 'cluster'),
        score_nb_of_hits = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('blastn', {}).get('score_nb_of_hits', 5)
    run:
        from camel.app.tools.pipelines.genedetection.blasthitfiltering import BlastHitFiltering
        hit_filtering = BlastHitFiltering(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(hit_filtering, input)
        step = Step(str(rule), hit_filtering, Camel.get_instance(), Path(str(params.running_dir)), wildcards)

        # Update parameters
        hit_filtering.update_parameters(
            min_percent_identity=str(params.min_percent_identity),
            min_coverage=str(params.min_coverage),
            filtering_method=str(params.filtering_method),
            score_nb_of_hits=str(params.score_nb_of_hits)
        )
        # Run tool
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(hit_filtering, output)

        # Add the informs from the filtering to the existing ones with the blastn command
        informs = SnakemakeUtils.load_object(Path(input.INFORMS_blastn))
        for key, value in hit_filtering.informs.items():
            if key.startswith('_'):
                continue
            informs[key] = value
        SnakemakeUtils.dump_object(informs, Path(output.INFORMS))

rule gene_detection_blast_text_alignment_generation:
    """
    Generates alignments in the text format.
    """
    input:
        ASN = rules.gene_detection_blast_blastn.output.ASN
    output:
        TXT = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'alignment_generation' / 'txt.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'alignment_generation'
    run:
        from camel.app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(blast_formatter, input)
        step = Step(str(rule), blast_formatter, Camel.get_instance(), Path(str(params.running_dir)), wildcards)
        blast_formatter.update_parameters(output_format='0', num_alignments=20000)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blast_formatter, output)

rule gene_detection_blast_text_alignment_extraction:
    """
    Extracts a text alignment for the selected hits and attaches them to the hit objects.
    """
    input:
        TXT = rules.gene_detection_blast_text_alignment_generation.output.TXT,
        VAL_Hits = rules.gene_detection_blast_hit_filtering.output.VAL_Hits
    output:
        VAL_Hits = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_HITS_METHOD).format(method='blast', db='{db}')
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'alignment_generation'
    run:
        from camel.app.tools.pipelines.genedetection.alignmentextractor import AlignmentExtractor
        alignment_extractor = AlignmentExtractor(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(alignment_extractor, input)
        step = Step(str(rule), alignment_extractor, Camel.get_instance(), Path(str(params.running_dir)), wildcards)
        step.run_step()
        hits_with_alignment = []
        for io_value, alignment in zip(
                SnakemakeUtils.load_object(Path(input.VAL_Hits)), alignment_extractor.tool_outputs['TXT']):
            io_value.value.alignment_path = Path(alignment.path)
            hits_with_alignment.append(io_value)
        SnakemakeUtils.dump_object(hits_with_alignment, Path(output.VAL_Hits))
