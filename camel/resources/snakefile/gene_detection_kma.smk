from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
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
        fasta_path = Path(SnakemakeUtils.load_object(Path(input.FASTA))[0].path)
        with open(fasta_path.parent / 'db_metadata.txt') as handle:
            metadata = json.load(handle)
        dir_kma = fasta_path.parent / 'kma'
        if not dir_kma.exists():
            raise FileNotFoundError(f"KMA database not found: {dir_kma}")
        kma_path = fasta_path.parent / 'kma' / fasta_path.stem
        SnakemakeUtils.dump_object([ToolIOValue(kma_path)], Path(output.DB))

rule gene_detection_kma:
    """
    Runs KMA on a database with the gene detection.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
        DB = rules.gene_detection_kma_get_db.output.DB
    output:
        TSV = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'kma' / 'tsv-kma.io',
        INFORMS = Path(config['working_dir']) / 'gene_detection' / '{db}' / 'kma' / 'informs_pre.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'kma',
        input_type = config.get('input_type', 'illumina'),
        ont = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('kma', {}).get('ont'),
        apm = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('kma', {}).get('apm'),
        cge = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('kma', {}).get('cge'),
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.kma.kma import KMA
        kma = KMA(Camel.get_instance())
        SnakemakeUtils.add_pickle_input(kma, 'DB', Path(input.DB))
        key_reads = 'PE' if params.input_type == 'illumina' else 'SE'
        fq_input_dict = SnakePipelineUtils.extracts_fq_input(
            Path(input.IO), key_pe='FASTQ_PE', key_se='FASTQ_SE', read_type=key_reads)
        kma.add_input_files(fq_input_dict)
        step = Step(str(rule), kma, Camel.get_instance(), Path(str(params.dir_)))
        if params.input_type == 'ont':
            kma.update_parameters(bc_nano=None, basecalls='0.7')
        if params.ont:
            kma.update_parameters(ont=None)
        if params.cge:
            kma.update_parameters(cge=None)
        if params.apm is not None:
            kma.update_parameters(apm=params.apm)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(kma, output)

rule gene_detection_kma_add_parameters:
    """
    Adds the filtering parameters to the informs (so that they are added to the HTML output report).
    """
    input:
        INFORMS = rules.gene_detection_kma.output.INFORMS
    output:
        INFORMS = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS_METHOD).format(method='kma', db='{db}')
    params:
        min_identity = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('kma',{}).get('min_percent_identity', 90),
        min_coverage = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('kma',{}).get('min_coverage', 60)
    run:
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        informs['Min. percent identity'] = f'{float(str(params.min_identity)):.0f}%'
        informs['Min. coverage'] = f'{float(str(params.min_coverage)):.0f}%'
        SnakemakeUtils.dump_object(informs, Path(output.INFORMS))

rule gene_detection_kma_hit_extraction:
    """
    Extracts and filters the hits detected by KMA.
    """
    input:
        TSV = rules.gene_detection_kma.output.TSV
    output:
        VAL_hits = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_HITS_METHOD).format(method='kma', db='{db}')
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'gene_detection' / wildcards.db / 'kma',
        min_identity = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('kma', {}).get('min_percent_identity', 90),
        min_coverage = lambda wildcards: config['gene_detection'][wildcards.db].get('params', {}).get('kma', {}).get('min_coverage', 60)
    run:
        from camel.app.tools.kma.kmagenedetectionhitextractor import KMAGeneDetectionHitExtractor
        extractor = KMAGeneDetectionHitExtractor(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(extractor, input)
        step = Step(str(rule), extractor, Camel.get_instance(), Path(str(params.dir_)))
        extractor.update_parameters(
            min_percent_identity=float(str(params.min_identity)),
            min_percent_coverage=float(str(params.min_coverage))
        )
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(extractor, output)
