import shutil
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
        CRAI = Path(config['final_output_dir']) / f'{config["sample"]}.cram.crai',

        VCF = Path(config['final_output_dir']) / f'{config["sample"]}.gVCF.gz',


rule move_output:
    input:
        CRAM = Path(config['working_dir']) / bam_to_cram.OUTPUT_BAMTOCRAM_CRAM,
        CRAM_checksum = Path(config['working_dir']) / bam_to_cram.OUTPUT_BAMTOCRAM_CRAM_checksum,
        CRAI = Path(config['working_dir']) / bam_to_cram.OUTPUT_BAMTOCRAM_CRAI,
        VCF = Path(config['working_dir']) / variant_calling.OUTPUT_gVCF,
        VCF_index = Path(config['working_dir']) / variant_calling.OUTPUT_gVCF_index,
        QC_done = Path(config['final_output_dir']) / 'qc' / 'qc_done.txt'
    output:
        CRAM = Path(config['final_output_dir']) / f'{config["sample"]}.cram',
        CRAM_checksum = Path(config['final_output_dir']) / f'{config["sample"]}.cram.md5',
        CRAI = Path(config['final_output_dir']) / f'{config["sample"]}.cram.crai',
        VCF = Path(config['final_output_dir']) / f'{config["sample"]}.gVCF.gz',
        VCF_index = Path(config['final_output_dir']) / f'{config["sample"]}.gVCF.gz.tbi',
    run:
        SnakemakeUtils.load_object(Path(input.CRAM))[0].path.link_to(output.CRAM)
        SnakemakeUtils.load_object(Path(input.CRAI))[0].path.link_to(output.CRAI)
        Path(input.CRAM_checksum).link_to(output.CRAM_checksum)
        SnakemakeUtils.load_object(Path(input.VCF))[0].path.link_to(output.VCF)
        Path(input.VCF_index).link_to(output.VCF_index)

        if not config['debug']:
            shutil.rmtree(Path(config['working_dir']) / 'alignment')
            shutil.rmtree(Path(config['working_dir']) / 'bamtocram')
            shutil.rmtree(Path(config['working_dir']) / 'input')
            shutil.rmtree(Path(config['working_dir']) / 'ref_input')
            shutil.rmtree(Path(config['working_dir']) / 'variant_calling')
            shutil.rmtree(Path(config['working_dir']) / 'qc')

rule prepare_references_io:
    """
    Prepare reference genome IO files
    """
    input:
        fasta_genome = Path(config['references']['ref_fasta']),
        dict_genome = Path(config['references']['ref_dict']),
        dbsnp = Path(config['references']['dbsnp_vcf']),
        known_indels = [Path(f) for f in config['references']['known_indels_sites_vcfs']],
        calling_intervals = Path(config['references']["calling_interval_list"]),
        contamination_sites_ud = Path(config['references']["contamination_sites_ud"]),
        contamination_sites_bed = Path(config['references']["contamination_sites_bed"]),
        contamination_sites_mu = Path(config['references']["contamination_sites_mu"]),
        coverage_interval_list = Path(config['references']["coverage_interval_list"]),
        evaluation_interval_list = Path(config['references']["evaluation_interval_list"]),
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
        SnakemakeUtils.dump_object([ToolIOValue(Path(input.fasta_genome))], Path(output.FASTA_GENOME))
        SnakemakeUtils.dump_object([ToolIOFile(Path(input.fasta_genome))], Path(output.FASTA_GENOME_FILE))
        SnakemakeUtils.dump_object([ToolIOFile(Path(input.dict_genome))], Path(output.DICT_GENOME))
        SnakemakeUtils.dump_object([ToolIOFile(Path(input.dbsnp))], Path(output.DBSNP))
        SnakemakeUtils.dump_object([ToolIOFile(Path(f)) for f in input.known_indels], Path(output.KNOWN_INDELS))
        SnakemakeUtils.dump_object([ToolIOFile(Path(input.calling_intervals))], Path(output.CALLING_INTERVALS))
        SnakemakeUtils.dump_object([ToolIOFile(Path(input.contamination_sites_bed)), ToolIOFile(Path(input.contamination_sites_mu)), ToolIOFile(Path(input.contamination_sites_ud))], Path(output.CONTAMINATION_SITES_UD))
        SnakemakeUtils.dump_object([ToolIOFile(Path(input.coverage_interval_list))], Path(output.COVERAGE_INTERVALS))
        SnakemakeUtils.dump_object([ToolIOFile(Path(input.evaluation_interval_list))], Path(output.EVALUATION_INTERVALS))

rule prepare_interval_pickles:
    input:
        interval_files = expand(Path(config['intervals_location']) / 'interval_{i}.intervals', i = config['intervals'])
    output:
        pickled_interval_files = expand(Path(config['working_dir']) / 'input' / 'interval_files' / 'interval_{i}.intervals.io', i = config['intervals'])
    params:
        output_dir = Path(config['working_dir']) / 'input' / 'interval_files'
    run:
        for interval_file in input.interval_files:
            SnakemakeUtils.dump_object([ToolIOFile(Path(interval_file))], Path(params.output_dir) / f"{Path(interval_file).name}.io")


rule move_qc:
    input:
        TXT = expand(Path(config['working_dir']) / "qc" / "quality_yield" / "{fastq}.unmapped.quality_yield_metrics.io", fastq = config["input_basenames"]),

        base_distribution_by_cycle = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.base_distribution_by_cycle.pdf", fastq = config["input_basenames"]),
        base_distribution_by_cycle_metrics = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.base_distribution_by_cycle_metrics", fastq = config["input_basenames"]),
        insert_size_histogram_pdf = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.insert_size_histogram.pdf", fastq = config["input_basenames"]),
        insert_size_metrics = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.insert_size_metrics", fastq = config["input_basenames"]),
        quality_by_cycle_pdf = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.quality_by_cycle.pdf", fastq = config["input_basenames"]),
        quality_by_cycle_metrics = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.quality_by_cycle_metrics", fastq = config["input_basenames"]),
        quality_distribution_pdf = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.quality_distribution.pdf", fastq = config["input_basenames"]),
        quality_distribution_metrics = expand(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.quality_distribution_metrics", fastq = config["input_basenames"]),

        alignment_summary_metrics = Path(config['working_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.alignment_summary_metrics",
        gc_bias_detail_metrics = Path(config['working_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.gc_bias.detail_metrics",
        gc_bias_pdf = Path(config['working_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.gc_bias.pdf",
        gc_bias_summary_metrics = Path(config['working_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.gc_bias.summary_metrics",

        mark_duplicates_metrics = Path(config['working_dir']) / 'qc' / 'mark_duplicates' / "duplicate_metrics.txt.io",

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

        TXT_metrics_CRAM = Path(config['working_dir']) / "qc" / "bamtocram" / "cram_validation_report.io",

        TXT_metrics_varCalling = Path(config['working_dir']) / "qc" / "variant_calling_metrics" / f"{config['sample']}.variant_calling_metrics.io",

        QC_summary = Path(config['working_dir']) / "qc" / "QC_summary.txt"

    output:
        TXT = expand(Path(config['final_output_dir']) / "qc" / "quality_yield" / "{fastq}.unmapped.quality_yield_metrics.txt", fastq = config["input_basenames"]),

        ## picard_unsorted_RG_quality
        base_distribution_by_cycle = expand(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.base_distribution_by_cycle.pdf", fastq = config["input_basenames"]),
        base_distribution_by_cycle_metrics = expand(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.base_distribution_by_cycle_metrics", fastq = config["input_basenames"]),
        insert_size_histogram_pdf = expand(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.insert_size_histogram.pdf", fastq = config["input_basenames"]),
        insert_size_metrics = expand(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.insert_size_metrics", fastq = config["input_basenames"]),
        quality_by_cycle_pdf = expand(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.quality_by_cycle.pdf", fastq = config["input_basenames"]),
        quality_by_cycle_metrics = expand(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.quality_by_cycle_metrics", fastq = config["input_basenames"]),
        quality_distribution_pdf = expand(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.quality_distribution.pdf", fastq = config["input_basenames"]),
        quality_distribution_metrics = expand(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / "{fastq}.unsorted_readgroup.quality_distribution_metrics", fastq = config["input_basenames"]),       
        
        alignment_summary_metrics = Path(config['final_output_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.alignment_summary_metrics",
        gc_bias_detail_metrics = Path(config['final_output_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.gc_bias.detail_metrics",
        gc_bias_pdf = Path(config['final_output_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.gc_bias.pdf",
        gc_bias_summary_metrics = Path(config['final_output_dir']) / "qc" / "RG_quality" / f"{config['sample']}.readgroup.gc_bias.summary_metrics",

        mark_duplicates_metrics = Path(config['final_output_dir']) / 'qc' / 'mark_duplicates' / f"{config['sample']}.duplicate_metrics.txt",

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
        
        TXT_metrics_CRAM = Path(config['final_output_dir']) / "qc" / "bamtocram" / f"{config['sample']}.cram_validation_report.txt",
        
        TXT_metrics_varCalling = Path(config['final_output_dir']) / "qc" / "variant_calling_metrics" / f"{config['sample']}.variant_calling_metrics.txt",

        QC_summary = Path(config['final_output_dir']) / "qc" / "QC_summary.txt",

        QC_done = Path(config['final_output_dir']) / 'qc' / 'qc_done.txt'
    run:
        for quality_yield_io in input.TXT:
            quality_yield_file = SnakemakeUtils.load_object(Path(quality_yield_io))[0].path
            file = Path(quality_yield_file).parts[-1]
            out = Path(config['final_output_dir']) / "qc" / "quality_yield" / f'{file}.txt'
            quality_yield_file.link_to(out)

        [Path(f).link_to(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / Path(f).name) for f in input.base_distribution_by_cycle]
        [Path(f).link_to(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / Path(f).name) for f in input.base_distribution_by_cycle_metrics]
        [Path(f).link_to(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / Path(f).name) for f in input.insert_size_histogram_pdf]
        [Path(f).link_to(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / Path(f).name) for f in input.insert_size_metrics]
        [Path(f).link_to(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / Path(f).name) for f in input.quality_by_cycle_pdf]
        [Path(f).link_to(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / Path(f).name) for f in input.quality_by_cycle_metrics]
        [Path(f).link_to(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / Path(f).name) for f in input.quality_distribution_metrics]
        [Path(f).link_to(Path(config['final_output_dir']) / "qc" / "unsorted_RG_quality" / Path(f).name) for f in input.quality_distribution_pdf]

        Path(input.alignment_summary_metrics).link_to(output.alignment_summary_metrics)
        Path(input.gc_bias_detail_metrics).link_to(output.gc_bias_detail_metrics)
        Path(input.gc_bias_pdf).link_to(output.gc_bias_pdf)
        Path(input.gc_bias_summary_metrics).link_to(output.gc_bias_summary_metrics)
        Path(input.alignment_summary_metrics_agg).link_to(output.alignment_summary_metrics_agg)
        Path(input.bait_bias_detail_metrics).link_to(output.bait_bias_detail_metrics)
        Path(input.bait_bias_summary_metrics).link_to(output.bait_bias_summary_metrics)
        SnakemakeUtils.load_object(Path(input.mark_duplicates_metrics))[0].path.link_to(output.mark_duplicates_metrics)
        Path(input.gc_bias_detail_metrics_agg).link_to(output.gc_bias_detail_metrics_agg)
        Path(input.gc_bias_pdf_agg).link_to(output.gc_bias_pdf_agg)
        Path(input.gc_bias_summary_metrics_agg).link_to(output.gc_bias_summary_metrics_agg)
        Path(input.insert_size_histogram_pdf_agg).link_to(output.insert_size_histogram_pdf_agg)
        Path(input.insert_size_metrics_agg).link_to(output.insert_size_metrics_agg)
        Path(input.pre_adapter_detail_metrics).link_to(output.pre_adapter_detail_metrics)
        Path(input.pre_adapter_summary_metrics).link_to(output.pre_adapter_summary_metrics)
        Path(input.quality_distribution_pdf_agg).link_to(output.quality_distribution_pdf_agg)
        Path(input.quality_distribution_metrics_agg).link_to(output.quality_distribution_metrics_agg)
        Path(input.error_summary_metrics).link_to(output.error_summary_metrics)
        Path(input.TXT_metrics_checksum).link_to(output.TXT_metrics_checksum)
        Path(input.TXT_metrics_WGS).link_to(output.TXT_metrics_WGS)
        Path(input.TXT_metrics_rawWGS).link_to(output.TXT_metrics_rawWGS)
        SnakemakeUtils.load_object(Path(input.TXT_metrics_CRAM))[0].path.link_to(output.TXT_metrics_CRAM)
        SnakemakeUtils.load_object(Path(input.TXT_metrics_varCalling))[0].path.link_to(output.TXT_metrics_varCalling)
        Path(input.QC_summary).link_to(output.QC_summary)

        subprocess.run(f"touch {output.QC_done}", shell = True, executable='/bin/bash')