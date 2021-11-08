from pathlib import Path

from camel.app.camel import Camel
from camel.app.components.pipelines import pipeutils
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

camel = Camel.get_instance()

rule bwa_aln_to_bam:
    input:
        FASTQ = Path(config['working_dir']) / 'input' / "{input_basename}.fastq.gz.io",
        FASTA_GENOME = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value.io",
    output:
        BAM = Path(config['working_dir']) / "alignment" / "mapping" / "{input_basename}.aligned.bam.io",
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / 'alignment' / "mapping" / wildcards.input_basename,
        output_file = lambda wildcards: f'{wildcards.input_basename}.aligned.unsorted.bam'
    threads: config["params_smk"]["threads_bwa"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_bwa"]
    run:
        from camel.app.tools.bwa.bwamap import BWAMap
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        Path(params.working_dir).mkdir(exist_ok=True)

        bwa_mem = BWAMap(camel)
        sam_to_bam = SamtoolsView(camel)

        SnakemakeUtils.add_pickle_input(bwa_mem, 'FASTQ_PE', Path(input.FASTQ))
        SnakemakeUtils.add_pickle_input(bwa_mem, 'INDEX_GENOME_PREFIX', Path(input.FASTA_GENOME))

        bwa_mem.update_parameters(
            threads = threads,
            **config['rule_params']['alignment'][rule]
        )

        sam_to_bam.update_parameters(
            output_filename = params.output_file,
            threads = threads
        )

        pipeutils.run_as_pipe([bwa_mem, sam_to_bam], params.working_dir)

        SnakemakeUtils.dump_tool_output(sam_to_bam, "BAM", Path(output.BAM))

rule picard_add_readgroups:
    input:
        BAM = rules.bwa_aln_to_bam.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'alignment' / 'add_readgroups' / '{input_basename}.aligned_rgadded.bam.io'
    params:
        working_dir = Path(config['working_dir']) / 'alignment' / 'add_readgroups',
        output_file = lambda wildcards: f'{wildcards.input_basename}.aligned_rgadded.bam',
        RG_id = lambda wildcards: f'{wildcards.input_basename}.rg'
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.addorreplacereadgroups import AddOrReplaceReadGroups

        add_rg = AddOrReplaceReadGroups(camel)
        SnakemakeUtils.add_pickle_inputs(add_rg, input)
        step = Step(rule, add_rg, camel, params.working_dir)
        add_rg.update_parameters(
            output = params.output_file,
            RG_id = params.RG_id,
            RG_sample_name = config['sample'],
            **config['rule_params']['alignment'][rule],
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(add_rg, "BAM", Path(output.BAM))

rule picard_mark_duplicates_sort:
    """
    Aggregate aligned+merged flowcell BAM files and mark duplicates
    We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
    to avoid having to spend time just merging BAM files.
    """
    input:
        BAM = expand(rules.picard_add_readgroups.output.BAM, input_basename = config['input_basenames'])
    output:
        BAM = Path(config['working_dir']) / 'alignment' / 'mark_duplicates' / f'{config["sample"]}.aligned.sorted.duplicates_marked.bam.io',
        metrics = Path(config['working_dir']) / 'qc' / 'mark_duplicates' / f"{config['sample']}.duplicate_metrics.txt.io"
    params:
        working_dir = Path(config['working_dir']) / 'alignment' / 'mark_duplicates',
        output_file = f'{config["sample"]}.aligned.unsorted.duplicates_marked.bam',
        metrics_output_file = Path(config['working_dir']) / 'qc' / 'mark_duplicates' / f"{config['sample']}.duplicate_metrics.txt"
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_mark_duplicates_sort"]
    run:
        from camel.app.tools.picard.markduplicates import MarkDuplicates
        from camel.app.tools.picard.sortsam import SortSam

        mark_duplicates = MarkDuplicates(camel)
        sort_sam = SortSam(camel)

        mark_duplicates.add_input_files({"BAM": [SnakemakeUtils.load_object(Path(path))[0] for path in input.BAM]})

        mark_duplicates.update_parameters(
            output = params.output_file,
            metrics_output = params.metrics_output_file,
            **config['rule_params']['alignment']['picard_mark_duplicates'],
        )
        sort_sam.update_parameters(
            output = "aligned.duplicate_marked.sorted.bam",
            **config['rule_params']['alignment']['picard_sort_sam'],
        )
        mark_duplicates.update_java_options(f"-mx{resources.mem_mb}M -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dpicard.useLegacyParser=false")
        sort_sam.update_java_options(f"-mx{resources.mem_mb}M -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dpicard.useLegacyParser=false")

        pipeutils.run_as_pipe([mark_duplicates, sort_sam], params.working_dir)

        SnakemakeUtils.dump_tool_output(mark_duplicates, "METRICS", Path(output.metrics))
        SnakemakeUtils.dump_tool_output(sort_sam, "BAM", Path(output.BAM))

rule picard_set_tags:
    input:
        BAM = rules.picard_mark_duplicates_sort.output.BAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io"
    output:
        BAM = Path(config['working_dir']) / 'alignment' / 'set_tags' / 'aligned_rgadded_dup-removed_settags.bam.io'
    params:
        working_dir = Path(config['working_dir']) / 'alignment' / 'set_tags',
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.setnmmdanduqtags import SetNmMdAndUqTags

        set_tags = SetNmMdAndUqTags(camel)
        SnakemakeUtils.add_pickle_inputs(set_tags, input)
        step = Step(rule, set_tags, camel, params.working_dir)
        set_tags.update_parameters(
            output = "aligned_dupmarked_sorted_rgadded_settags.bam",
            **config['rule_params']['alignment'][rule]
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(set_tags, "BAM", Path(output.BAM))

rule gatk4_baserecalibrator:
    input:
        BAM = rules.picard_set_tags.output.BAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
        VCF_KNOWN_SNPS = Path(config['working_dir']) / "ref_input" / "dbsnp_vcf.io",
        VCF_KNOWN_INDELS = Path(config['working_dir']) / "ref_input" / "known_indels_vcf.io",
        TXT_intervals = Path(config['working_dir']) / 'input' / 'interval_files' / "interval_{i}.intervals.io"
    output:
        TXT_RecalibrationTable = Path(config['working_dir']) / "alignment" / "bqsr" / "{i}_recal_data.csv.io"
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / "alignment" / "bqsr",
        output_file = lambda wildcards: f'{wildcards.i}_recal_data.csv'
    threads: config["params_smk"]["threads_bqsr"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_bqsr"]
    run:
        from camel.app.tools.gatk4.gatk4baserecalibrator import GATK4BaseRecalibrator

        bqsr = GATK4BaseRecalibrator(camel)
        SnakemakeUtils.add_pickle_input(bqsr, "BAM", Path(input.BAM))
        SnakemakeUtils.add_pickle_input(bqsr, "FASTA_REF", Path(input.FASTA_REF))
        SnakemakeUtils.add_pickle_input(bqsr, "VCF_KNOWN_SNPS", Path(input.VCF_KNOWN_SNPS))
        SnakemakeUtils.add_pickle_input(bqsr, "VCF_KNOWN_INDELS", Path(input.VCF_KNOWN_INDELS))
        SnakemakeUtils.add_pickle_input(bqsr, "TXT_intervals", Path(input.TXT_intervals))

        step = Step(rule, bqsr, camel, params.working_dir)
        bqsr.update_parameters(
            **config['rule_params']['alignment'][rule],
            output = params.output_file
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(bqsr, "TXT_RecalibrationTable", Path(output.TXT_RecalibrationTable))

rule gatk4_gather_bqsr_reports:
    input:
        bqsr_report_interval = expand(rules.gatk4_baserecalibrator.output.TXT_RecalibrationTable, i = config["intervals"])
    output:
        bqsr_report_gathered = Path(config['working_dir']) / "alignment" / "gather_bqsr_reports" / "recal_data.csv.io",
    params:
        working_dir = Path(config['working_dir']) / "alignment" / "gather_bqsr_reports",
    threads: config["params_smk"]["threads_bqsr"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_bqsr"]
    run:
        from camel.app.tools.gatk4.gatk4gatherbqsrreports import GATK4GatherBQSRReports

        Path(params.working_dir).mkdir(exist_ok=True)

        gather_bqsr = GATK4GatherBQSRReports(camel)
        step = Step(rule, gather_bqsr, camel, params.working_dir)
        gather_bqsr.add_input_files({"TXT_intervals": [SnakemakeUtils.load_object(Path(path))[0] for path in input.bqsr_report_interval]})
        gather_bqsr.update_parameters(
            output = "recal_data.csv"
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(gather_bqsr, 'TXT_RecalibrationTable', Path(output.bqsr_report_gathered))

rule gatk4_apply_bqsr:
    input:
        BAM = rules.gatk4_baserecalibrator.input.BAM, # Ensure it gets same input as baserecalibrator
        BQSR = rules.gatk4_gather_bqsr_reports.output.bqsr_report_gathered,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
        TXT_intervals = Path(config['working_dir']) / 'input' / 'interval_files' / "interval_{i}.intervals.io"
    output:
        BAM = Path(config['working_dir']) / "alignment" / "apply_bqsr" / "{i}_sorted.bam.io"
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / "alignment" / "apply_bqsr",
        output_file = lambda wildcards: f'{wildcards.i}_sorted.bam'
    threads: config["params_smk"]["threads_apply_bqsr"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_bqsr"]
    run:
        from camel.app.tools.gatk4.gatk4applybqsr import GATK4ApplyBQSR

        Path(params.working_dir).mkdir(exist_ok=True)

        apply_bqsr = GATK4ApplyBQSR(camel)
        SnakemakeUtils.add_pickle_inputs(apply_bqsr, input)
        step = Step(rule, apply_bqsr, camel, params.working_dir)
        apply_bqsr.update_parameters(
            **config['rule_params']['alignment'][rule],
            output = params.output_file
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(apply_bqsr, "BAM", Path(output.BAM))

rule picard_gather_sorted_bam:
    input:
        bqsr_BAM_interval = expand(rules.gatk4_apply_bqsr.output.BAM, i = config["intervals"])
    output:
        BAM = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam" / "bqsr_gathered_sorted.bam.io",
    params:
        working_dir = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam"
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.gatherbamfiles import GatherBamFiles

        Path(params.working_dir).mkdir(exist_ok=True)

        gather_bam = GatherBamFiles(camel)
        step = Step(rule, gather_bam, camel, params.working_dir)
        gather_bam.add_input_files({"BAMs": [SnakemakeUtils.load_object(Path(path))[0] for path in input.bqsr_BAM_interval]})
        gather_bam.update_parameters(
            **config['rule_params']['alignment'][rule],
            output = "gathered_sorted.bam"
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(gather_bam, 'BAM', Path(output.BAM))
