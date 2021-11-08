from pathlib import Path

from bidict import bidict

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.broadwgs.snakefile import alignment, bam_to_cram, variant_calling

camel = Camel.get_instance()

rule picard_quality_yield:
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_INTERMEDIATE_BAM,
    output:
        TXT = Path(config['working_dir']) / "qc" / "quality_yield" / "{input_basename}.unmapped.quality_yield_metrics.io"
    params:
        working_dir = Path(config['working_dir']) / "qc" / "quality_yield",
        output_file = lambda wildcards: f"{wildcards.input_basename}.unmapped.quality_yield_metrics"
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectqualityyieldmetrics import CollectQualityYieldMetrics

        quality_yield = CollectQualityYieldMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(quality_yield, input)
        step = Step(rule, quality_yield, camel, params.working_dir)
        quality_yield.update_parameters(**config['rule_params']['qc'][rule], output = params.output_file)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quality_yield, output)

rule picard_unsorted_RG_quality:
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_INTERMEDIATE_BAM,
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
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectmultiplemetrics import CollectMultipleMetrics

        unsorted_RG_quality = CollectMultipleMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(unsorted_RG_quality, input)
        step = Step(rule, unsorted_RG_quality, camel, params.working_dir)
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
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectmultiplemetrics import CollectMultipleMetrics

        RG_quality = CollectMultipleMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(RG_quality, input)
        step = Step(rule, RG_quality, camel, params.working_dir)
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
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectmultiplemetrics import CollectMultipleMetrics

        agg_metrics = CollectMultipleMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(agg_metrics, input)
        step = Step(rule, agg_metrics, camel, params.working_dir)
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
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.calculatereadgroupchecksum import CalculateReadGroupChecksum

        get_checksum = CalculateReadGroupChecksum(camel)
        SnakemakeUtils.add_pickle_input(get_checksum, "BAM", Path(input.BAM))
        step = Step(rule, get_checksum, camel, params.working_dir)
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
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectwgsmetrics import CollectWgsMetrics

        Path(params.working_dir).mkdir(exist_ok=True)

        wgs_metrics = CollectWgsMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(wgs_metrics, input)
        step = Step(rule, wgs_metrics, camel, params.working_dir)
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
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectrawwgsmetrics import CollectRawWgsMetrics

        raw_wgs_metrics = CollectRawWgsMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(raw_wgs_metrics, input)
        step = Step(rule, raw_wgs_metrics, camel, params.working_dir)
        raw_wgs_metrics.update_parameters(
            output = params.output_file,
            **config['rule_params']['qc'][rule]
        )
        step.run_step()


rule picard_validate_cram:
    input:
        CRAM = Path(config['working_dir']) / bam_to_cram.OUTPUT_BAMTOCRAM_CRAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        TXT_metrics = Path(config['working_dir']) / "qc" / "bamtocram" / "cram_validation_report.io",
    params:
        working_dir = Path(config['working_dir']) / "qc" / "bamtocram"
    threads: config["params_smk"]["threads_cram"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_cram"]
    run:
        from camel.app.tools.picard.validatesamfile import ValidateSamFile

        Path(params.working_dir).mkdir(exist_ok=True)

        val_cram = ValidateSamFile(camel)
        SnakemakeUtils.add_pickle_input(val_cram, "BAM", Path(input.CRAM))
        SnakemakeUtils.add_pickle_input(val_cram, "FASTA_REF", Path(input.FASTA_REF))
        step = Step(rule, val_cram, camel, params.working_dir)
        val_cram.update_parameters(
            **config['rule_params']['bam_to_cram'][rule]
        )
        val_cram.update_java_options("-mx100G -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dpicard.useLegacyParser=false")
        step.run_step()
        SnakemakeUtils.dump_tool_output(val_cram, 'TXT_report', Path(output.TXT_metrics))

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
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectvariantcallingmetrics import CollectVariantCallingMetrics

        Path(params.working_dir).mkdir(exist_ok=True)

        vcf_metrics = CollectVariantCallingMetrics(camel)
        SnakemakeUtils.add_pickle_inputs(vcf_metrics, input)
        step = Step(rule, vcf_metrics, camel, params.working_dir)
        vcf_metrics.update_parameters(
            output_prefix = config["sample"],
            **config['rule_params']['qc'][rule]
        )
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(vcf_metrics, output)

rule generate_qc_summary:
    input:
        mark_duplicates_metrics = alignment.OUTPUT_MARK_DUPLICATES_METRICS,
        alignment_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f'{config["sample"]}.agg.alignment_summary_metrics',
        insert_size_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f'{config["sample"]}.agg.insert_size_metrics',
        wgs_metrics = rules.picard_wgs_metrics.output.TXT_metrics,
        variant_calling_metrics = rules.picard_variant_calling_metrics.output.TXT_report
    output:
        QC_summary = Path(config['working_dir']) / "qc" / "QC_summary.txt"
    run:
        alignment_metrics = {}
        variant_calling_metrics = {}

        # Mark Duplicates metrics
        with open(input.mark_duplicates_metrics) as handle:
            lines = handle.readlines()
            metrics_class_index = lines.index("## METRICS CLASS\tpicard.sam.DuplicationMetrics\n")

        header = lines[metrics_class_index + 1].split("\t")
        values = lines[metrics_class_index + 2].split("\t")

        mark_duplicates = dict(zip(header, values))
        alignment_metrics['percent_duplication'] = mark_duplicates['PERCENT_DUPLICATION']

        # Alignment summary metrics
        with open("/home/chdevogelaere/qc/aggregation_metrics/GC088815.agg.alignment_summary_metrics") as handle:
            lines = handle.readlines()
            metrics_class_index = lines.index("## METRICS CLASS\tpicard.analysis.AlignmentSummaryMetrics\n")

        header = lines[metrics_class_index + 1].split("\t")
        values = lines[metrics_class_index + 4].split("\t")

        agg_alignment_summary = dict(zip(header, values))
        alignment_metrics['percent_chimeras'] = agg_alignment_summary['PCT_CHIMERAS']

        # Insert size metrics
        with open("/home/chdevogelaere/qc/aggregation_metrics/GC088815.agg.insert_size_metrics") as handle:
            lines = handle.readlines()
            metrics_class_index = lines.index("## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\n")

        header = lines[metrics_class_index + 1].split("\t")
        values = lines[metrics_class_index + 2].split("\t")

        agg_insert_size = dict(zip(header, values))
        alignment_metrics['median_insert_size'] = agg_insert_size['MEDIAN_INSERT_SIZE']

        # WGS metrics
        with open("/home/chdevogelaere/qc/wgs_metrics/GC088815.wgs.metrics.txt") as handle:
            lines = handle.readlines()
            metrics_class_index = lines.index("## METRICS CLASS\tpicard.analysis.WgsMetrics\n")
            histogram_start_index = lines.index("## HISTOGRAM\tjava.lang.Integer\n")

        header = lines[metrics_class_index + 1].split("\t")
        values = lines[metrics_class_index + 2].split("\t")

        wgs_metrics = dict(zip(header, values))
        mean_coverage = float(wgs_metrics['MEAN_COVERAGE'])
        genome_territory = int(wgs_metrics['GENOME_TERRITORY'])
        lowest_20 = float(int(genome_territory) * 0.2)

        coverage_histo = lines[histogram_start_index+3:-1]
        coverage_histo = [l.rstrip().split("\t") for l in coverage_histo]

        cov_values = [int(l[1]) for l in coverage_histo]
        histo_cumsum = [sum(cov_values[:i+1]) for i, value in enumerate(cov_values)]
        histo_cumsum_bidict = bidict([(i, n) for i, n in enumerate(histo_cumsum)])

        tmp = histo_cumsum + [int(lowest_20)]
        tmp.sort()
        find_index = tmp.index(int(lowest_20))
        val_before = histo_cumsum[find_index - 1]
        val_after = histo_cumsum[find_index]

        cov_20 = (((val_after - lowest_20) / (val_after - val_before)) + find_index - 1)

        fold80 = mean_coverage / cov_20
        alignment_metrics['evenness_coverage'] = fold80

        # Variant calling metrics
        with open("/home/chdevogelaere/qc/variant_calling_metrics/GC088815.variant_calling_metrics.txt") as handle:
            lines = handle.readlines()
            metrics_class_index = lines.index("## METRICS CLASS\tpicard.vcf.CollectVariantCallingMetrics$VariantCallingSummaryMetrics\n")

        header = lines[metrics_class_index + 1].split("\t")
        values = lines[metrics_class_index + 2].split("\t")

        variant_calling_summary = dict(zip(header, values))
        variant_calling_metrics["TITV_DBSNP"] = variant_calling_summary['DBSNP_TITV']
        variant_calling_metrics["TITV_NOVEL"] = variant_calling_summary['NOVEL_TITV']
