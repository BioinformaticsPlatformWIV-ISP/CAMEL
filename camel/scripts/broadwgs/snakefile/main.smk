import os
import subprocess

from pathlib import Path

from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
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
        CRAM = Path(config['final_output_dir']) / f'{config["sample"]}.cram',
        CRAM_checksum = Path(config['final_output_dir']) / f'{config["sample"]}.cram.md5',
        CRAI = Path(config['final_output_dir']) / f'{config["sample"]}.crai',
        CRAM_metrics = Path(config['final_output_dir']) / f'{config["sample"]}.cram.metrics',

        VCF = Path(config['final_output_dir']) / f'{config["sample"]}.gVCF',

        QC_done = Path(config['final_output_dir']) / 'qc' / 'qc_done.txt'


rule move_output:
    input:
        CRAM = Path(config['working_dir']) / bam_to_cram.OUTPUT_BAMTOCRAM_CRAM,
        CRAM_checksum = Path(config['working_dir']) / bam_to_cram.OUTPUT_BAMTOCRAM_CRAM_checksum,
        CRAI = Path(config['working_dir']) / bam_to_cram.OUTPUT_BAMTOCRAM_CRAI,
        CRAM_metrics = Path(config['working_dir']) / bam_to_cram.OUTPUT_BAMTOCRAM_CRAM_metrics,
        VCF = Path(config['working_dir']) / variant_calling.OUTPUT_gVCF,
    output:
        CRAM = Path(config['final_output_dir']) / f'{config["sample"]}.cram',
        CRAM_checksum = Path(config['final_output_dir']) / f'{config["sample"]}.cram.md5',
        CRAI = Path(config['final_output_dir']) / f'{config["sample"]}.crai',
        CRAM_metrics = Path(config['final_output_dir']) / f'{config["sample"]}.cram.metrics',
        VCF = Path(config['final_output_dir']) / f'{config["sample"]}.gVCF',
    run:
        os.link(SnakemakeUtils.load_object(input.CRAM)[0].path, output.CRAM)
        os.link(SnakemakeUtils.load_object(input.CRAI)[0].path, output.CRAI)
        os.link(input.CRAM_checksum, output.CRAM_checksum)
        os.link(SnakemakeUtils.load_object(input.CRAM_metrics)[0].path, output.CRAM_metrics)
        os.link(SnakemakeUtils.load_object(input.VCF)[0].path, output.VCF)


rule prepare_references_io:
    """
    Prepare reference genome IO files
    """
    input:
        fasta_genome = config['references']['ref_fasta'],
        dict_genome = config['references']['ref_dict'],
        dbsnp = config['references']['dbsnp_vcf'],
        known_indels = config['references']['known_indels_sites_vcfs'],
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
        KNOWN_INDELS = Path(config['working_dir']) / "ref_input" / "known_indels_vcf.io",
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
        io_knownINDELs = [ToolIOFile(f) for f in input.known_indels]
        io_calling_intervals = [ToolIOFile(input.calling_intervals)]
        io_contamination = [ToolIOFile(input.contamination_sites_bed), ToolIOFile(input.contamination_sites_mu), ToolIOFile(input.contamination_sites_ud)]
        io_coverage_intervals = [ToolIOFile(input.coverage_interval_list)]
        io_evaluation_intervals = [ToolIOFile(input.evaluation_interval_list)]

        # Dump objects
        SnakemakeUtils.dump_object(io_fasta_genome, str(output.FASTA_GENOME))
        SnakemakeUtils.dump_object(io_fasta_genome_file, str(output.FASTA_GENOME_FILE))
        SnakemakeUtils.dump_object(io_dict_genome, str(output.DICT_GENOME))
        SnakemakeUtils.dump_object(io_dbSNP, str(output.DBSNP))
        SnakemakeUtils.dump_object(io_knownINDELs, str(output.KNOWN_INDELS))
        SnakemakeUtils.dump_object(io_calling_intervals, str(output.CALLING_INTERVALS))
        SnakemakeUtils.dump_object(io_contamination, str(output.CONTAMINATION_SITES_UD))
        SnakemakeUtils.dump_object(io_coverage_intervals, str(output.COVERAGE_INTERVALS))
        SnakemakeUtils.dump_object(io_evaluation_intervals, str(output.EVALUATION_INTERVALS))

rule move_qc:
    input:
        alignment_summary_metrics = Path(config['working_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.alignment_summary_metrics",
        gc_bias_detail_metrics = Path(config['working_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.gc_bias.detail_metrics",
        gc_bias_pdf = Path(config['working_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.gc_bias.pdf",
        gc_bias_summary_metrics = Path(config['working_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.gc_bias.summary_metrics",

        alignment_summary_metrics_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.alignment_summary_metrics",
        bait_bias_detail_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.bait_bias_detail_metrics",
        bait_bias_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.bait_bias_summary_metrics",
        gc_bias_detail_metrics_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.gc_bias.detail_metrics",
        gc_bias_pdf_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.gc_bias.pdf",
        gc_bias_summary_metrics_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.gc_bias.summary_metrics",
        insert_size_histogram_pdf_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.insert_size_histogram.pdf",
        insert_size_metrics_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.insert_size_metrics",
        pre_adapter_detail_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.pre_adapter_detail_metrics",
        pre_adapter_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.pre_adapter_summary_metrics",
        quality_distribution_pdf_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.quality_distribution.pdf",
        quality_distribution_metrics_agg = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.quality_distribution_metrics",
        error_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.error_summary_metrics",

        TXT_metrics_checksum = Path(config['working_dir']) / "qc" / "RG_checksum" / f"{config['sample']}.bam.read_group_md5",

        TXT_metrics_WGS = Path(config['working_dir']) / "qc" / "wgs_metrics" / f"{config['sample']}.wgs.metrics.txt",
        TXT_metrics_rawWGS = Path(config['working_dir']) / "qc" / "wgs_metrics" / f"{config['sample']}.raw.wgs.metrics.txt",

        TXT_metrics_validateGVCF = Path(config['working_dir']) / "qc" / "validate_gvcf" / f"{config['sample']}.validate_vcf.txt",

        TXT_metrics_varCalling = Path(config['working_dir']) / "qc" / "variant_calling_metrics" / f"{config['sample']}.variant_calling_metrics.io",

        ## TODO uBAM metrics
        # quality_yield = expand(Path(config['working_dir']) / "qc" / "ubam_quality_yield" / "{fastq}.unmapped.quality_yield_metrics.io", fastq = config["input_basenames"]),
        # base_distribution_by_cycle = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.base_distribution_by_cycle.pdf", fastq = config["input_basenames"]),
        # base_distribution_by_cycle_metrics = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.base_distribution_by_cycle_metrics", fastq = config["input_basenames"]),
        # insert_size_histogram_pdf = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.insert_size_histogram.pdf", fastq = config["input_basenames"]),
        # insert_size_metrics = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.insert_size_metrics", fastq = config["input_basenames"]),
        # quality_by_cycle_pdf = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.quality_by_cycle.pdf", fastq = config["input_basenames"]),
        # quality_by_cycle_metrics = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.quality_by_cycle_metrics", fastq = config["input_basenames"]),
        # quality_distribution_pdf = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.quality_distribution.pdf", fastq = config["input_basenames"]),
        # quality_distribution_metrics = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.quality_distribution_metrics", fastq = config["input_basenames"]),
    output:
        alignment_summary_metrics = Path(config['final_output_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.alignment_summary_metrics",
        gc_bias_detail_metrics = Path(config['final_output_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.gc_bias.detail_metrics",
        gc_bias_pdf = Path(config['final_output_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.gc_bias.pdf",
        gc_bias_summary_metrics = Path(config['final_output_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.gc_bias.summary_metrics",

        alignment_summary_metrics_agg = Path(config['final_output_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.alignment_summary_metrics",
        bait_bias_detail_metrics = Path(config['final_output_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.bait_bias_detail_metrics",
        bait_bias_summary_metrics = Path(config['final_output_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.bait_bias_summary_metrics",
        gc_bias_detail_metrics_agg = Path(config['final_output_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.gc_bias.detail_metrics",
        gc_bias_pdf_agg = Path(config['final_output_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.gc_bias.pdf",
        gc_bias_summary_metrics_agg = Path(config['final_output_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.gc_bias.summary_metrics",
        insert_size_histogram_pdf_agg = Path(config['final_output_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.insert_size_histogram.pdf",
        insert_size_metrics_agg = Path(config['final_output_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.insert_size_metrics",
        pre_adapter_detail_metrics = Path(config['final_output_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.pre_adapter_detail_metrics",
        pre_adapter_summary_metrics = Path(config['final_output_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.pre_adapter_summary_metrics",
        quality_distribution_pdf_agg = Path(config['final_output_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.quality_distribution.pdf",
        quality_distribution_metrics_agg = Path(config['final_output_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.quality_distribution_metrics",
        error_summary_metrics = Path(config['final_output_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.error_summary_metrics",

        TXT_metrics_checksum = Path(config['final_output_dir']) / "qc" / "RG_checksum" / f"{config['sample']}.bam.read_group_md5",

        TXT_metrics_WGS = Path(config['final_output_dir']) / "qc" / "wgs_metrics" / f"{config['sample']}.wgs.metrics.txt",
        TXT_metrics_rawWGS = Path(config['final_output_dir']) / "qc" / "wgs_metrics" / f"{config['sample']}.raw.wgs.metrics.txt",

        TXT_metrics_validateGVCF = Path(config['final_output_dir']) / "qc" / "validate_gvcf" / f"{config['sample']}.validate_vcf.txt",

        TXT_metrics_varCalling = Path(config['final_output_dir']) / "qc" / "variant_calling_metrics" / f"{config['sample']}.variant_calling_metrics.txt",

        QC_done = Path(config['final_output_dir']) / 'qc' / 'qc_done.txt'
    run:
        os.link(input.alignment_summary_metrics, output.alignment_summary_metrics)
        os.link(input.gc_bias_detail_metrics, output.gc_bias_detail_metrics)
        os.link(input.gc_bias_pdf, output.gc_bias_pdf)
        os.link(input.gc_bias_summary_metrics, output.gc_bias_summary_metrics)
        os.link(input.alignment_summary_metrics_agg, output.alignment_summary_metrics_agg)
        os.link(input.bait_bias_detail_metrics, output.bait_bias_detail_metrics)
        os.link(input.bait_bias_summary_metrics, output.bait_bias_summary_metrics)
        os.link(input.gc_bias_detail_metrics_agg, output.gc_bias_detail_metrics_agg )
        os.link(input.gc_bias_pdf_agg, output.gc_bias_pdf_agg)
        os.link(input.gc_bias_summary_metrics_agg, output.gc_bias_summary_metrics_agg)
        os.link(input.insert_size_histogram_pdf_agg, output.insert_size_histogram_pdf_agg)
        os.link(input.insert_size_metrics_agg, output.insert_size_metrics_agg)
        os.link(input.pre_adapter_detail_metrics, output.pre_adapter_detail_metrics)
        os.link(input.pre_adapter_summary_metrics, output.pre_adapter_summary_metrics)
        os.link(input.quality_distribution_pdf_agg, output.quality_distribution_pdf_agg )
        os.link(input.quality_distribution_metrics_agg, output.quality_distribution_metrics_agg)
        os.link(input.error_summary_metrics, output.error_summary_metrics)
        os.link(input.TXT_metrics_checksum, output.TXT_metrics_checksum)
        os.link(input.TXT_metrics_WGS, output.TXT_metrics_WGS)
        os.link(input.TXT_metrics_rawWGS, output.TXT_metrics_rawWGS)
        os.link(input.TXT_metrics_validateGVCF, output.TXT_metrics_validateGVCF)
        os.link(SnakemakeUtils.load_object(input.TXT_metrics_varCalling)[0].path, output.TXT_metrics_varCalling)

        subprocess.run(f"touch {output.QC_done}", shell = True, executable='/bin/bash')