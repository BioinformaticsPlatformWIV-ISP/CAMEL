from pathlib import Path

from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import gene_detection

rule gene_detection_kma_get_db:
    """
    Retrieves the database for running KMA.
    """
    input:
        FASTA = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'db_manager' / 'fasta-clust.io'
    output:
        DB = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'kma' / 'db.io'
    run:
        import json
        from camel.app.io.tooliovalue import ToolIOValue
        fasta_path = Path(SnakemakeUtils.load_object(input.FASTA)[0].path)
        with open(fasta_path.parent / 'db_metadata.txt') as handle:
            metadata = json.load(handle)
        dir_kma = fasta_path.parent / 'kma'
        if not dir_kma.exists():
            raise FileNotFoundError(f"KMA database not found: {dir_kma}")
        kma_path = fasta_path.parent / 'kma' / fasta_path.stem
        SnakemakeUtils.dump_object([ToolIOValue(kma_path)], output.DB)

rule gene_detection_kma:
    """
    Runs KMA on a database with the gene detection.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
        DB = rules.gene_detection_kma_get_db.output.DB
    output:
        TSV = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'kma' / 'tsv-kma.io',
        INFORMS = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS_METHOD).format(method='kma', db='{db}')
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'kma'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.kma.kma import KMA
        kma = KMA(camel)
        SnakemakeUtils.add_pickle_input(kma, 'DB', input.DB)
        kma.add_input_files(SnakePipelineUtils.extracts_fq_input(input.IO))
        step = Step(rule, kma, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(kma, output)

rule gene_detection_kma_hit_extraction:
    """
    Extracts and filters the hits detected by KMA.
    """
    input:
        TSV = rules.gene_detection_kma.output.TSV
    output:
        VAL_hits = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_HITS_METHOD).format(method='kma', db='{db}')
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'kma',
        min_identity = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('kma', {}).get('min_percent_identity', 90),
        min_coverage = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('kma', {}).get('min_coverage', 60)
    run:
        from camel.app.tools.kma.kmagenedetectionhitextractor import KMAGeneDetectionHitExtractor
        extractor = KMAGeneDetectionHitExtractor(camel)
        SnakemakeUtils.add_pickle_inputs(extractor, input)
        step = Step(rule, extractor, camel, params.running_dir, config)
        extractor.update_parameters(
            min_percent_identity=float(params.min_identity),
            min_percent_coverage=float(params.min_coverage)
        )
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(extractor, output)
