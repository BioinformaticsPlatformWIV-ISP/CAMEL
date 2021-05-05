from pathlib import Path

from camel.scripts.broadwgs.snakefile import alignment, bam_to_cram, variant_calling, qc

#######################
# Included Snakefiles #
#######################
include: alignment.SNAKEFILE_ALIGNMENT
include: bam_to_cram.SNAKEFILE_BAMTOCRAM
include: variant_calling.SNAKEFILE_VARIANTCALLING
include: qc.SNAKEFILE_QC

#########
# Rules #
#########

rule all:
    input:
        BAM = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam" / (config["sample"] + ".bam.io"),

        CRAM = Path(config['working_dir']) / "bamtocram" / "convert" / (config["sample"] + ".cram.io"),
        CRAM_checksum = Path(config['working_dir']) / "bamtocram" / (config["sample"] + ".cram.md5"),
        CRAI = Path(config['working_dir']) / "bamtocram" / "index" / (config["sample"] + ".sample.crai.io"),
        CRAM_metrics = Path(config['working_dir']) / "bamtocram" / "metrics" / (config["sample"] + "_cram_validation_report.io"),

        VCF = Path(config['working_dir']) / "variant_calling" / "merge_vcf" / (config["sample"] + ".g.vcf.gz.io"),

        #QC_done = Path(config['working_dir']) / 'qc' / 'qc_done.txt'


rule prepare_references_io:
    """
    Prepare reference genome IO files for snakemake to use: generate io files for
    """
    input:
        fasta_genome = config['references']['ref_fasta'],
        dict_genome = config['references']['ref_dict'],
        dbsnp = config['references']['dbsnp_vcf'],
        #known_indels
        calling_intervals = config['references']["calling_interval_list"],
        contamination_sites_ud = config['references']["contamination_sites_ud"],
        contamination_sites_bed = config['references']["contamination_sites_bed"],
        contamination_sites_mu = config['references']["contamination_sites_mu"],
        coverage_interval_list = config['references']["coverage_interval_list"],
        evaluation_interval_list = config['references']["evaluation_interval_list"],
    output:
        FASTA_GENOME = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value.io",
        FASTA_GENOME_FILE = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
        DICT_GENOME = Path(config['working_dir']) / "ref_input" / "dictionary_genome_human.io",
        DBSNP = Path(config['working_dir']) / "ref_input" / "dbsnp_vcf.io",
        CALLING_INTERVALS = Path(config['working_dir']) / "ref_input" / "calling_intervals.io",
        CONTAMINATION_SITES_UD = Path(config['working_dir']) / "ref_input" / "contamination_sites.io",
        COVERAGE_INTERVALS = Path(config['working_dir']) / "ref_input" / "coverage_interval_list.io",
        EVALUATION_INTERVALS = Path(config['working_dir']) / "ref_input" / "evaluation_interval_list.io",
    run:
        # Make objects
        io_fasta_genome = [ToolIOValue(input.fasta_genome)]
        io_fasta_genome_file = [ToolIOFile(input.fasta_genome)]
        io_dict_genome = [ToolIOFile(input.dict_genome)]
        io_dbSNP = [ToolIOFile(input.dbsnp)]
        io_calling_intervals = [ToolIOFile(input.calling_intervals)]
        io_contamination = [ToolIOFile(input.contamination_sites_bed), ToolIOFile(input.contamination_sites_mu), ToolIOFile(input.contamination_sites_ud)]
        io_coverage_intervals = [ToolIOFile(input.coverage_interval_list)]
        io_evaluation_intervals = [ToolIOFile(input.evaluation_interval_list)]

        # Dump objects
        SnakemakeUtils.dump_object(io_fasta_genome, str(output.FASTA_GENOME))
        SnakemakeUtils.dump_object(io_fasta_genome_file, str(output.FASTA_GENOME_FILE))
        SnakemakeUtils.dump_object(io_dict_genome, str(output.DICT_GENOME))
        SnakemakeUtils.dump_object(io_dbSNP, str(output.DBSNP))
        SnakemakeUtils.dump_object(io_calling_intervals, str(output.CALLING_INTERVALS))
        SnakemakeUtils.dump_object(io_contamination, str(output.CONTAMINATION_SITES_UD))
        SnakemakeUtils.dump_object(io_coverage_intervals, str(output.COVERAGE_INTERVALS))
        SnakemakeUtils.dump_object(io_evaluation_intervals, str(output.EVALUATION_INTERVALS))

rule run_qc:
    input:
        TXT = expand(Path(config['working_dir']) / "qc" / "ubam_quality_yield" / config["sample"] / "{ubam}.unmapped.quality_yield_metrics.io", ubam = config["ubams"]),

        ## picard_unsorted_RG_quality
        base_distribution_by_cycle = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.base_distribution_by_cycle.pdf", ubam = config["ubams"]),
        base_distribution_by_cycle_metrics = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.base_distribution_by_cycle_metrics", ubam = config["ubams"]),
        insert_size_histogram_pdf = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.insert_size_histogram.pdf", ubam = config["ubams"]),
        insert_size_metrics = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.insert_size_metrics", ubam = config["ubams"]),
        quality_by_cycle_pdf = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.quality_by_cycle.pdf", ubam = config["ubams"]),
        quality_by_cycle_metrics = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.quality_by_cycle_metrics", ubam = config["ubams"]),
        quality_distribution_pdf = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.quality_distribution.pdf", ubam = config["ubams"]),
        quality_distribution_metrics = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{ubam}.unsorted_readgroup.quality_distribution_metrics", ubam = config["ubams"]),

        ## picard_RG_quality
        alignment_summary_metrics = Path(config['working_dir']) / "qc" / "RG_quality" / (config["sample"] + ".readgroup.alignment_summary_metrics"),
        gc_bias_detail_metrics = Path(config['working_dir']) / "qc" / "RG_quality" / (config["sample"] + ".readgroup.gc_bias.detail_metrics"),
        gc_bias_pdf = Path(config['working_dir']) / "qc" / "RG_quality" / (config["sample"] + ".readgroup.gc_bias.pdf"),
        gc_bias_summary_metrics = Path(config['working_dir']) / "qc" / "RG_quality" / (config["sample"] + ".readgroup.gc_bias.summary_metrics"),

        ## picard_aggregation_metrics
        alignment_summary_metrics_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.alignment_summary_metrics"),
        bait_bias_detail_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.bait_bias_detail_metrics"),
        bait_bias_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.bait_bias_summary_metrics"),
        gc_bias_detail_metrics_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.gc_bias.detail_metrics"),
        gc_bias_pdf_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.gc_bias.pdf"),
        gc_bias_summary_metrics_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.gc_bias.summary_metrics"),
        insert_size_histogram_pdf_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.insert_size_histogram.pdf"),
        insert_size_metrics_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.insert_size_metrics"),
        pre_adapter_detail_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.pre_adapter_detail_metrics"),
        pre_adapter_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.pre_adapter_summary_metrics"),
        quality_distribution_pdf_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.quality_distribution.pdf"),
        quality_distribution_metrics_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.quality_distribution_metrics"),
        error_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / (config["sample"] + ".agg.error_summary_metrics"),

        TXT_metrics_checksum = Path(config['working_dir']) / "qc" / "RG_checksum" / (config["sample"] + ".bam.read_group_md5.io"),
        TXT_metrics_WGS = Path(config['working_dir']) / "qc" / "wgs_metrics" / (config["sample"] + ".wgs.metrics.io"),
        TXT_metrics_rawWGS = Path(config['working_dir']) / "qc" / "wgs_metrics" / (config["sample"] + ".raw.wgs.metrics.io"),
        TXT_metrics_validateGVCF = Path(config['working_dir']) / "qc" / "validate_gvcf" / (config["sample"] + ".validate_vcf.io"),
        TXT_metrics_varCalling = Path(config['working_dir']) / "qc" / "variant_calling_metrics" / (config["sample"] + ".variant_calling_metrics.io"),
    output:
        done = Path(config['working_dir']) / 'qc' / 'qc_done.txt'
    shell:
        "touch {output.done}"

