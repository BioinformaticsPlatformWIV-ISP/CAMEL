from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.influenzapipeline.snakefile import alignment
from camel.scripts.influenzapipeline.snakefile import genometyping_blastn


camel = Camel.get_instance()


rule run_alignment:
    """
    Aligns reads to a given reference genome using BWA or Bowtie2. 
    
    Input is FASTQ_PE in case the rule is run during the iterative consensus sequence calling. 
    In that case, FASTQ_PE will be present in the config file.
    """
    input:
        FASTQ_PE = config.get('FASTQ_PE', []),
        IO = Path(config['working_dir']) / 'fq_dict.io' if 'FASTQ_PE' not in config else [],
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_INDEX_GENOME_PREFIX if 'index_genome_prefix' not in config else config['index_genome_prefix']
    output:
        SAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_SAM,
        INFORMS = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'alignment'
    threads: 6
    run:
        from camel.app.tools.bwa.bwamap import BWAMap
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        mapper_class = BWAMap if config['aligner'] == 'bwa' else Bowtie2Map
        mapper = mapper_class(camel)

        if input.IO:
            fq_dict = SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE', key_se='FASTQ_SE')
            input_files = fq_dict['FASTQ_PE']
            mapper.add_input_files({'FASTQ_PE': input_files})
            SnakemakeUtils.add_pickle_input(mapper, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        else:
            SnakemakeUtils.add_pickle_inputs(mapper, input, keys=['FASTQ_PE', 'INDEX_GENOME_PREFIX'])

        step = Step(str(rule), mapper, camel, params.running_dir)
        if config['aligner'] == 'bowtie2':
            if 'mapping' in config['rule_parameters']:
                mapper.update_parameters(**{'sensitive': False, 'end_to_end': False})
                mapper.update_parameters(**{k:True for k in config['rule_parameters']['mapping']['bowtie2_mode']})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(mapper, output)

rule alignment_sort_sam:
    input:
        SAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_SAM
    output:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM
    params:
        running_dir = Path(config['working_dir']) / 'alignment'
    threads: 6
    run:
        from camel.app.tools.picard.sortsam import SortSam

        sortsam = SortSam(camel)
        SnakemakeUtils.add_pickle_inputs(sortsam, input)
        step = Step(str(rule), sortsam, camel, params.running_dir)
        sortsam.update_parameters(**{'create_index': 'true', 'sort_order': 'coordinate'})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sortsam, output)

rule alignment_collect_metrics:
    input:
        FASTA_REF = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_FASTA_REF,
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM
    output:
        INFORMS = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_PICARD_METRICS,
        TXT_AlignmentSummary = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_ALIGNMENTSUMMARY,
        TXT_InsertSize = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_INSERTSIZE,
        TXT_QualityDistribution = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_MAPQUALITYDISTRIBUTION,
        TXT_QualityDistributionFigure = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_MAPQUALITYDISTRIBUTION_PDF,
        TXT_GcBias = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_GCBIAS,
        TXT_GcBiasSummary = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_GCBIAS_SUMMARY,
        TXT_GcBiasFigure = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_GCBIAS_FIGURE
    params:
        running_dir = Path(config['working_dir']) / 'alignment'
    threads: 6
    run:
        from camel.app.tools.picard.collectmultiplemetrics import CollectMultipleMetrics

        cmm = CollectMultipleMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(cmm, input)
        cmm_params = {'output_prefix': 'readmap_qc', 'reset_metrics': 'null',
                      'metrics_CollectAlignmentSummaryMetrics': 'CollectAlignmentSummaryMetrics',
                      'metrics_CollectInsertSizeMetrics': 'CollectInsertSizeMetrics',
                      'metrics_QualityScoreDistribution': 'QualityScoreDistribution',
                      'metrics_CollectGcBiasMetrics': 'CollectGcBiasMetrics'}
        step = Step(str(rule), cmm, camel, params.running_dir)
        cmm.update_parameters(**cmm_params)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(cmm, output)

rule alignment_samtools_depth:
    input:
        BAM =Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM
    output:
        TSV = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_SAMTOOLS_DEPTH,
        INFORMS = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_SAMTOOLS_DEPTH_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'alignment'
    threads: 6
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth

        sd = SamtoolsDepth(camel)
        SnakemakeUtils.add_pickle_inputs(sd, input)
        step = Step(str(rule), sd, camel, params.running_dir)
        sd.update_parameters(**{'maximum_coverage_depth': 1000000})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sd, output)

rule alignment_samtools_depth_analyzer:
    input:
        TXT = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_SAMTOOLS_DEPTH,
         FASTA_REF = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_FASTA_REF
    output:
        INFORMS = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_SAMTOOLS_DEPTH_ANALYZER_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'alignment'
    threads: 6
    run:
        from camel.app.tools.samtools.samtoolsdepthstatsanalyzer import SamtoolsDepthStatsAnalyzer

        sda = SamtoolsDepthStatsAnalyzer(camel)
        SnakemakeUtils.add_pickle_inputs(sda, input)
        step = Step(str(rule), sda, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sda, output)

rule alignment_report:
    input:
        INFORMS_alignment = rules.run_alignment.output.INFORMS,
        INFORMS_picardmetrics = rules.alignment_collect_metrics.output.INFORMS,
        INFORMS_samtoolsdepth = rules.alignment_samtools_depth_analyzer.output.INFORMS,
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
        PDF_GC = rules.alignment_collect_metrics.output.TXT_GcBiasFigure,
        PDF_MQC = rules.alignment_collect_metrics.output.TXT_QualityDistributionFigure
    output:
        VAL_HTML = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'alignment' / 'report'
    run:
        from camel.app.tools.pipelines.alignment.reporteralignment import ReporterAlignment

        reporter = ReporterAlignment(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, camel, params.running_dir)
        if 'genome_segments' in config['species_info']:
            reporter.update_parameters(**{'genome_segments': config['species_info']['genome_segments']})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)
