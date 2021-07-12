from pathlib import Path
import subprocess

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

from camel.scripts.broadwgs.snakefile import alignment, variant_calling

camel = Camel.get_instance()

rule picard_quality_yield:
    input:
        uBAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM_UNMAPPED,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        TXT = Path(config['working_dir']) / "qc" / "ubam_quality_yield" / "{input_basename}.unmapped.quality_yield_metrics.io"
    params:
        working_dir = Path(config['working_dir']) / "qc" / "ubam_quality_yield",
        output_file = lambda wildcards: f"{wildcards.input_basename}.unmapped.quality_yield_metrics"
    run:
        from camel.app.tools.picard.collectqualityyieldmetrics import CollectQualityYieldMetrics

        quality_yield = CollectQualityYieldMetrics(camel)
        SnakemakeUtils.add_pickle_input(quality_yield, "BAM", input.uBAM)
        SnakemakeUtils.add_pickle_input(quality_yield, "FASTA_REF", input.FASTA_REF)
        step = Step(rule, quality_yield, camel, params.working_dir, config)
        quality_yield.update_parameters(**config['rule_params']['qc'][rule])
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quality_yield, output)

rule picard_unsorted_RG_quality:
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM_UNMAPPED,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        multiext(str(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{input_basename}.unsorted_readgroup"),
                ".base_distribution_by_cycle.pdf",
                ".base_distribution_by_cycle_metrics",
                ".quality_by_cycle.pdf",
                ".quality_by_cycle_metrics",
                ".insert_size_histogram.pdf",
                ".insert_size_metrics",
                ".quality_distribution.pdf",
                ".quality_distribution_metrics"
                )

    params:
        working_dir = Path(config['working_dir']) / "qc" / "unsorted_RG_quality",
        output_prefix = lambda wildcards: f"{wildcards.input_basename}.unsorted_readgroup"
    run:
        from camel.app.tools.picard.collectmultiplemetrics import CollectMultipleMetrics

        unsorted_RG_quality = CollectMultipleMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(unsorted_RG_quality, input)
        step = Step(rule, unsorted_RG_quality, camel, params.working_dir, config)
        unsorted_RG_quality.update_parameters(
            output_prefix = params.output_prefix,
            **config['rule_params']['qc'][rule]
        )
        step.run_step()

rule picard_RG_quality:
    input:
        BAM = Path(config["working_dir"]) / alignment.OUTPUT_ALIGNMENT_BAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        multiext(str(Path(config['working_dir']) / "qc" / "RG_quality" / f'{config["sample"]}.readgroup'),
        ".alignment_summary_metrics",
        ".gc_bias.detail_metrics",
        ".gc_bias.pdf",
        ".gc_bias.summary_metrics")
    params:
        working_dir = Path(config['working_dir']) / "qc" / "RG_quality",
        output_prefix = f'{config["sample"]}.readgroup'
    run:
        from camel.app.tools.picard.collectmultiplemetrics import CollectMultipleMetrics

        RG_quality = CollectMultipleMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(RG_quality, input)
        step = Step(rule, RG_quality, camel, params.working_dir, config)
        RG_quality.update_parameters(
            output_prefix = params.output_prefix,
            **config['rule_params']['qc'][rule]
        )
        step.run_step()

# TODO maybe implement this later
# rule check_prevalidation:
#     input:
#         duplicate_metrics = Path(config['working_dir']) / alignment.OUTPUT_MARK_DUPLICATES_METRICS,
#         chimerism_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.alignment_summary_metrics",
#     output:
#         duplication_rate = "duplication_value.txt",
#         chimerism_rate = "chimerism_value.txt"
#     params:
#         working_dir = Path(config['working_dir']) / "qc" / "check_prevalidation",
#         max_duplication_in_reasonable_sample = 0.30,
#         max_chimerism_in_reasonable_sample = 0.15,
#         #Boolean is_outlier_data = duplication_rate > max_duplication_in_reasonable_sample || chimerism_rate > max_chimerism_in_reasonable_sample
#     run:
#         import subprocess
#
#         subprocess.run(f"grep -A 1 PERCENT_DUPLICATION {input.duplication_metrics} > {params.working_dir}.duplication.csv")
#         subprocess.run(f"grep -A 3 PCT_CHIMERAS {input.chimerism_metrics} | grep -v OF_PAIR > {params.working_dir}.chimerism.csv")
#
#         import csv
#         with open(f"{params.working_dir}.duplication.csv") as dupfile:
#             reader = csv.DictReader(dupfile, delimiter='\t')
#             for row in reader:
#                 with open(f"{params.working_dir}/{output.duplication_rate}","w") as file:
#                     file.write(row['PERCENT_DUPLICATION'])
#                     file.close()
#
#         with open(f"{params.working_dir}.chimerism.csv") as chimfile:
#             reader = csv.DictReader(chimfile, delimiter='\t')
#             for row in reader:
#                 with open(f"{params.working_dir}/{output.chimerism_rate}","w") as file:
#                     file.write(row['PCT_CHIMERAS'])
#                     file.close()


rule picard_aggregation_metrics:
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        multiext(str(Path(config['working_dir']) / "qc" / "aggregation_metrics" / f'{config["sample"]}.agg'),
        ".alignment_summary_metrics",
        ".bait_bias_detail_metrics",
        ".bait_bias_summary_metrics",
        ".gc_bias.detail_metrics",
        ".gc_bias.pdf",
        ".gc_bias.summary_metrics",
        ".insert_size_histogram.pdf",
        ".insert_size_metrics",
        ".pre_adapter_detail_metrics",
        ".pre_adapter_summary_metrics",
        ".quality_distribution.pdf",
        ".quality_distribution_metrics",
        ".error_summary_metrics")
    params:
          working_dir = Path(config['working_dir']) / "qc" / "aggregation_metrics",
          output_prefix = f'{config["sample"]}.agg',
    run:
        from camel.app.tools.picard.collectmultiplemetrics import CollectMultipleMetrics

        agg_metrics = CollectMultipleMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(agg_metrics, input)
        step = Step(rule, agg_metrics, camel, params.working_dir, config)
        agg_metrics.update_parameters(
            output_prefix = params.output_prefix,
            **config['rule_params']['qc'][rule]
        )
        step.run_step()

rule picard_RG_checksum:
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        TXT_checksum = Path(config['working_dir']) / "qc" / "RG_checksum" / f'{config["sample"]}.bam.read_group_md5',
    params:
        working_dir = Path(config['working_dir']) / "qc" / "RG_checksum",
        output_file = f'{config["sample"]}.bam.read_group_md5'
    run:
        from camel.app.tools.picard.calculatereadgroupchecksum import CalculateReadGroupChecksum

        get_checksum = CalculateReadGroupChecksum(camel)
        SnakemakeUtils.add_pickle_input(get_checksum, "BAM", input.BAM)
        step = Step(rule, get_checksum, camel, params.working_dir, config)
        get_checksum.update_parameters(output = params.output_file)
        step.run_step()


rule picard_wgs_metrics:
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
        COVERAGE_INTERVALS = Path(config['working_dir']) / "ref_input" / "coverage_interval_list.io"
    output:
        TXT_metrics = Path(config['working_dir']) / "qc" / "wgs_metrics" / f'{config["sample"]}.wgs.metrics.txt'
    params:
        working_dir = Path(config['working_dir']) / "qc" / "wgs_metrics",
        output_file = f'{config["sample"]}.wgs.metrics.txt'
    run:
        from camel.app.tools.picard.collectwgsmetrics import CollectWgsMetrics

        Path(params.working_dir).mkdir(exist_ok=True)

        wgs_metrics = CollectWgsMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(wgs_metrics, input)
        step = Step(rule, wgs_metrics, camel, params.working_dir, config)
        wgs_metrics.update_parameters(
            output = params.output_file,
            **config['rule_params']['qc'][rule]
        )
        step.run_step()

rule picard_raw_wgs_metrics:
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
        COVERAGE_INTERVALS = Path(config['working_dir']) / "ref_input" / "coverage_interval_list.io"
    output:
        TXT_metrics = Path(config['working_dir']) / "qc" / "wgs_metrics" / f'{config["sample"]}.raw.wgs.metrics.txt'
    params:
        working_dir = Path(config['working_dir']) / "qc" / "wgs_metrics",
        output_file = f'{config["sample"]}.raw.wgs.metrics.txt'
    run:
        from camel.app.tools.picard.collectrawwgsmetrics import CollectRawWgsMetrics

        raw_wgs_metrics = CollectRawWgsMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(raw_wgs_metrics, input)
        step = Step(rule, raw_wgs_metrics, camel, params.working_dir, config)
        raw_wgs_metrics.update_parameters(
            output = params.output_file,
            **config['rule_params']['qc'][rule]
        )
        step.run_step()

rule gatk4_validate_gvcf:
    input:
        VCF = Path(config['working_dir']) / variant_calling.OUTPUT_gVCF,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
        CALLING_INTERVALS = Path(config['working_dir']) / "ref_input" / "calling_intervals.io",
        VCF_dbsnp = Path(config['working_dir']) / "ref_input" / "dbsnp_vcf.io",
    output:
        TXT_metrics = Path(config['working_dir']) / "qc" / "validate_gvcf" / f'{config["sample"]}.validate_vcf.txt'
    params:
        working_dir = Path(config['working_dir']) / "qc" / "validate_gvcf",
        output_file = f'{config["sample"]}.validate_vcf.txt'
    run:
        from camel.app.tools.gatk4.gatk4validatevariants import GATK4ValidateVariants

        Path(params.working_dir).mkdir(exist_ok=True)

        validate_gvcf = GATK4ValidateVariants(camel)
        SnakemakeUtils.add_pickle_inputs(validate_gvcf, input)
        step = Step(rule, validate_gvcf, camel, params.working_dir, config)
        validate_gvcf.update_parameters(
            output = params.output_file,
            **config['rule_params']['qc'][rule]
        )
        step.run_step()


rule picard_variant_calling_metrics:
    input:
        VCF = Path(config['working_dir']) / variant_calling.OUTPUT_gVCF,
        VCF_dbsnp = Path(config['working_dir']) / "ref_input" / "dbsnp_vcf.io",
        DICT_GENOME = Path(config['working_dir']) / "ref_input" / "dictionary_genome_human.io",
        EVALUATION_INTERVALS = Path(config['working_dir']) / "ref_input" / "evaluation_interval_list.io",
    output:
        TXT_report = Path(config['working_dir']) / "qc" / "variant_calling_metrics" / f'{config["sample"]}.variant_calling_metrics.io',
    params:
        working_dir = Path(config['working_dir']) / "qc" / "variant_calling_metrics",
    run:
        from camel.app.tools.picard.collectvariantcallingmetrics import CollectVariantCallingMetrics

        Path(params.working_dir).mkdir(exist_ok=True)

        vcf_metrics = CollectVariantCallingMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(vcf_metrics, input)
        step = Step(rule, vcf_metrics, camel, params.working_dir, config)
        vcf_metrics.update_parameters(
            output_prefix = config["sample"],
            **config['rule_params']['qc'][rule]
        )
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(vcf_metrics, output)