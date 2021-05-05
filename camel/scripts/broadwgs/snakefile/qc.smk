from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue

camel = Camel()

rule picard_quality_yield:
    input:
        uBAM = Path(config['working_dir']) / "input" / config["sample"]  / "{ubam}.unmapped.bam.io"
    output:
        TXT = Path(config['working_dir']) / "qc" / "ubam_quality_yield" / config["sample"] / "{ubam}.unmapped.quality_yield_metrics.io"
    params:
        working_dir = Path(config['working_dir']) / "qc" / "ubam_quality_yield" / config["sample"]
    run:
        from camel.app.tools.picard.collectqualityyieldmetrics import CollectQualityYieldMetrics

        Path(params.working_dir).mkdir(exist_ok=True)

        quality_yield = CollectQualityYieldMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(quality_yield, input)
        step = Step(rule, quality_yield, camel, params.working_dir, config)
        quality_yield.update_parameters(original_qualities="true")
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quality_yield, output)

rule picard_unsorted_RG_quality:
    input:
        BAM = Path(config['working_dir']) / "alignment" / "merge_bam" / "{ubam}.aligned.unsorted.io",
    output:
        base_distribution_by_cycle = Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.base_distribution_by_cycle.pdf",
        base_distribution_by_cycle_metrics = Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.base_distribution_by_cycle_metrics",
        insert_size_histogram_pdf = Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.insert_size_histogram.pdf",
        insert_size_metrics = Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.insert_size_metrics",
        quality_by_cycle_pdf = Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.quality_by_cycle.pdf",
        quality_by_cycle_metrics = Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.quality_by_cycle_metrics",
        quality_distribution_pdf = Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.quality_distribution.pdf",
        quality_distribution_metrics = Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.quality_distribution_metrics"
    params:
        output_prefix = str(Path(config['working_dir']) / "qc" / "unsorted_RG_quality") + "{ubam}"
    shell:
        """
        run_picard.sh \
        CollectMultipleMetrics \
        INPUT={input.BAM} \
        OUTPUT={params.output_prefix} \
        ASSUME_SORTED=true \
        PROGRAM=null \
        PROGRAM=CollectBaseDistributionByCycle \
        PROGRAM=CollectInsertSizeMetrics \
        PROGRAM=MeanQualityByCycle \
        PROGRAM=QualityScoreDistribution \
        METRIC_ACCUMULATION_LEVEL=null \
        METRIC_ACCUMULATION_LEVEL=ALL_READS
        
        touch {params.output_prefix}.insert_size_metrics
        touch {params.output_prefix}.insert_size_histogram.pdf
        """

# TODO
# rule verify_bam_id:
#
#
#
# rule check_contamination:
#     input:
#         BAM_SORTED = Path(config['working_dir']) / "alignment" / "{sample}.aligned.duplicate_marked.sorted.bam.io",
#         FASTA_GENOME_FILE = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
#
#     output:
#         metrics = config['working_dir'] + "qc/{sample}.preBqsr.selfSM.io"
#         # output: read float from stout
#     params:
#         output_prefix = a,
#         contamination_underestimation_factor = 0.75
#     script:
#          "scripts/check_contamination.py"


rule picard_RG_quality:
    input:
        BAM = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam" / (config["sample"] + ".bam.io"),
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        alignment_summary_metrics = Path(config['working_dir']) / "qc" / "RG_quality" / (config["sample"]+ ".readgroup.alignment_summary_metrics"),
        gc_bias_detail_metrics = Path(config['working_dir']) / "qc" / "RG_quality" / (config["sample"] + ".readgroup.gc_bias.detail_metrics"),
        gc_bias_pdf = Path(config['working_dir']) / "qc" / "RG_quality" / (config["sample"] + ".readgroup.gc_bias.pdf"),
        gc_bias_summary_metrics = Path(config['working_dir']) / "qc" / "RG_quality" / (config["sample"] + ".readgroup.gc_bias.summary_metrics"),
    params:
        output_prefix = str(Path(config['working_dir']) / "qc" / "RG_quality") + (config["sample"] + ".readgroup")
    shell:
        """
        touch {params.output_prefix}.gc_bias.detail_metrics
        touch {params.output_prefix}.gc_bias.pdf
        touch {params.output_prefix}.gc_bias.summary_metrics
        
        run_picard.sh \
        CollectMultipleMetrics \
        INPUT={input.bam} \
        OUTPUT={params.output_prefix} \
        REFERENCE_SEQUENCE={input.ref_fasta} \
        ASSUME_SORTED=true \
        PROGRAM=null \
        PROGRAM=CollectAlignmentSummaryMetrics \
        PROGRAM=CollectGcBiasMetrics \
        METRIC_ACCUMULATION_LEVEL=null \
        METRIC_ACCUMULATION_LEVEL=READ_GROUP

        """


## TODO
#
# rule check_prevalidation:
#     input:
#         duplicate_metrics = config['working_dir'] + "qc/{sample}.duplicate_metrics",
#         chimerism_metrics = config['working_dir'] + "aggbamqc/{sample}.alignment_summary_metrics",
#     output:
#         duplication_rate = "duplication_value.txt", ##TODO
#         chimerism_rate = "chimerism_value.txt" ##TODO
#     params:
#         max_duplication_in_reasonable_sample = 0.30,
#         max_chimerism_in_reasonable_sample = 0.15,
#         #Boolean is_outlier_data = duplication_rate > max_duplication_in_reasonable_sample || chimerism_rate > max_chimerism_in_reasonable_sample
#     run:
#         import subprocess
#
#         subprocess.run(f"grep -A 1 PERCENT_DUPLICATION {input.duplication_metrics} > duplication.csv")
#         subprocess.run(f"grep -A 3 PCT_CHIMERAS {input.chimerism_metrics} | grep -v OF_PAIR > chimerism.csv")
#
#         import csv
#         with open({output.duplication_rate}) as dupfile:
#             reader = csv.DictReader(dupfile, delimiter='\t')
#             for row in reader:
#                 with open({"duplication.csv"},"w") as file:
#                     file.write(row['PERCENT_DUPLICATION'])
#                     file.close()
#
#         with open('chimerism.csv') as chimfile:
#             reader = csv.DictReader(chimfile, delimiter='\t')
#             for row in reader:
#                 with open({output.chimerism_rate},"w") as file:
#                     file.write(row['PCT_CHIMERAS'])
#                     file.close()


rule picard_aggregation_metrics:
    input:
        BAM = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam" / (config["sample"] + ".bam.io"),
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        alignment_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.alignment_summary_metrics"),
        bait_bias_detail_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.bait_bias_detail_metrics"),
        bait_bias_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.bait_bias_summary_metrics"),
        gc_bias_detail_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.gc_bias.detail_metrics"),
        gc_bias_pdf = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.gc_bias.pdf"),
        gc_bias_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.gc_bias.summary_metrics"),
        insert_size_histogram_pdf = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.insert_size_histogram.pdf"),
        insert_size_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.insert_size_metrics"),
        pre_adapter_detail_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.pre_adapter_detail_metrics"),
        pre_adapter_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.pre_adapter_summary_metrics"),
        quality_distribution_pdf = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.quality_distribution.pdf"),
        quality_distribution_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.quality_distribution_metrics"),
        error_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.error_summary_metrics"),
    params:
          output_prefix = str(Path(config['working_dir']) / "qc" / "aggregation_metrics") + config["sample"] + ".agg",
    shell:
        """
        touch {params.output_prefix}.gc_bias.detail_metrics \
              {params.output_prefix}.gc_bias.pdf \
              {params.output_prefix}.gc_bias.summary_metrics \
              {params.output_prefix}.insert_size_metrics \
              {params.output_prefix}.insert_size_histogram.pdf \
              {params.output_prefix}.quality_distribution.pdf

        run_picard.sh \
          CollectMultipleMetrics \
          INPUT={input.bam} \
          REFERENCE_SEQUENCE={input.ref_fasta} \
          OUTPUT={params.output_prefix} \
          ASSUME_SORTED=true \
          PROGRAM=null \
          PROGRAM=CollectAlignmentSummaryMetrics \
          PROGRAM=CollectInsertSizeMetrics \
          PROGRAM=CollectSequencingArtifactMetrics \
          PROGRAM=QualityScoreDistribution \
          PROGRAM=CollectGcBiasMetrics \
          METRIC_ACCUMULATION_LEVEL=null \
          METRIC_ACCUMULATION_LEVEL=SAMPLE \
          METRIC_ACCUMULATION_LEVEL=LIBRARY
        """

rule picard_RG_checksum:
    input:
        BAM = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam" / (config["sample"] + ".bam.io"),
    output:
        TXT_checksum = Path(config['working_dir']) / "qc" / "RG_checksum" / (config["sample"] + ".bam.read_group_md5.io"),
    params:
        working_dir = Path(config['working_dir']) / "qc" / "RG_checksum"
    run:
        from camel.app.tools.pîcard.calculatereadgroupchecksum import CalculateReadGroupChecksum

        Path(params.working_dir).mkdir(exist_ok=True)

        get_checksum = CalculateReadGroupChecksum(camel)
        SnakemakeUtils.add_pickle_input(get_checksum, "BAM", input.BAM)
        step = Step(rule, get_checksum, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(get_checksum, 'TXT_checksum', output.TXT_checksum)


rule picard_wgs_metrics:
    input:
        BAM = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam" / (config["sample"] + ".bam.io"),
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
        COVERAGE_INTERVALS = Path(config['working_dir']) / "ref_input" / "coverage_interval_list.io"
    output:
        TXT_metrics = Path(config['working_dir']) / "qc" / "wgs_metrics" / (config["sample"] + ".wgs.metrics.io")
    params:
        working_dir = Path(config['working_dir']) / "qc" / "wgs_metrics",
        read_length = config["parameters"]["read_length"]
    run:
        from camel.app.tools.picard.collectwgsmetrics import CollectWgsMetrics

        Path(params.working_dir).mkdir(exist_ok=True)

        wgs_metrics = CollectWgsMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(wgs_metrics, input)
        step = Step(rule, wgs_metrics, camel, params.working_dir, config)
        wgs_metrics.update_parameters(
            validation_stringency = "SILENT",
            include_bq_histogram = "true",
            use_fast_algorithm = "true",
            read_length = params.read_length
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(wgs_metrics, 'TXT_metrics', output.TXT_metrics)

"""
Collect raw WGS metrics
"""

rule picard_raw_wgs_metrics:
    input:
        BAM = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam" / (config["sample"] + ".bam.io"),
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
        COVERAGE_INTERVALS = Path(config['working_dir']) / "ref_input" / "coverage_interval_list.io"
    output:
        TXT_metrics = Path(config['working_dir']) / "qc" / "wgs_metrics" / (config["sample"] + ".raw.wgs.metrics.io")
    params:
        working_dir = Path(config['working_dir']) / "qc" / "wgs_metrics",
        read_length = config["parameters"]["read_length"]
    run:
        from camel.app.tools.picard.collectwgsmetrics import CollectRawWgsMetrics

        Path(params.working_dir).mkdir(exist_ok=True)

        raw_wgs_metrics = CollectWgsMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(raw_wgs_metrics, input)
        step = Step(rule, raw_wgs_metrics, camel, params.working_dir, config)
        raw_wgs_metrics.update_parameters(
            validation_stringency = "SILENT",
            include_bq_histogram = "true",
            use_fast_algorithm = "true",
            read_length = params.read_length
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(raw_wgs_metrics, 'TXT_metrics', output.TXT_metrics)

rule gatk4_validate_gvcf:
    input:
        VCF = Path(config['working_dir']) / "variant_calling" / "merge_vcf" / (config["sample"] + ".g.vcf.gz.io"),
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
        CALLING_INTERVALS = Path(config['working_dir']) / "ref_input" / "calling_intervals.io",
        VCF_dbsnp = Path(config['working_dir']) / "ref_input" / "dbsnp_vcf.io",
    output:
        TXT_metrics = Path(config['working_dir']) / "qc" / "validate_gvcf" / (config["sample"] + ".validate_vcf.io")
    params:
        working_dir = Path(config['working_dir']) / "qc" / "validate_gvcf"
    run:
        from camel.app.tools.gatk4.gatk4validatevariants import GATK4ValidateVariants

        Path(params.working_dir).mkdir(exist_ok=True)

        validate_gvcf = CollectWgsMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(validate_gvcf, input)
        step = Step(rule, validate_gvcf, camel, params.working_dir, config)
        validate_gvcf.update_parameters(
            gvcf = "true",
            val_type_to_exclude = "ALLELES"
        )
        step.run_step()
        #TODO no output -> informs?


rule picard_variant_calling_metrics:
    input:
        VCF = Path(config['working_dir']) / "variant_calling" / "merge_vcf" / (config["sample"] + ".g.vcf.gz.io"),
        VCF_dbsnp = Path(config['working_dir']) / "ref_input" / "dbsnp_vcf.io",
        DICT_GENOME = Path(config['working_dir']) / "ref_input" / "dictionary_genome_human.io",
        EVALUATION_INTERVALS = Path(config['working_dir']) / "ref_input" / "evaluation_interval_list.io",
    output:
        TXT_metrics = Path(config['working_dir']) / "qc" / "variant_calling_metrics" / (config["sample"] + ".variant_calling_metrics.io"),
    params:
        working_dir = Path(config['working_dir']) / "qc" / "variant_calling_metrics"
    run:
        from camel.app.tools.picard.collectvariantcallingmetrics import CollectVariantCallingMetrics

        Path(params.working_dir).mkdir(exist_ok=True)

        vcf_metrics = CollectVariantCallingMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(vcf_metrics, input)
        step = Step(rule, vcf_metrics, camel, params.working_dir, config)
        vcf_metrics.update_parameters(gvcf = true)
        step.run_step()
        SnakemakeUtils.dump_tool_output(vcf_metrics, 'TXT_metrics', output.TXT_metrics)