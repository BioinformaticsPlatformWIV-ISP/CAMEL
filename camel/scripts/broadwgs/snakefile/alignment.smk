from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.components.pipelines import pipeutils

camel = Camel.get_instance()

rule bwa_aln_to_bam:
    input:
        FASTQ = Path(config['working_dir']) / 'input' / "{input_basename}.fastq.gz.io",
        FASTA_GENOME = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value.io",
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io"
    output:
        BAM = Path(config['working_dir']) / "alignment" / "mapping" / "{input_basename}.aligned.bam.io",
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / 'alignment' / "mapping" / wildcards.input_basename,
        output_file = lambda wildcards: f'{wildcards.input_basename}.aligned.unsorted.bam'
    threads: config["params_smk"]["threads_bwa"]
    run:
        from camel.app.tools.bwa.bwamap import BWAMap
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        Path(params.working_dir).mkdir(exist_ok=True)

        bwa_mem = BWAMap(camel)
        sam_to_bam = SamtoolsView(camel)

        SnakemakeUtils.add_pickle_input(bwa_mem, 'FASTQ_PE', input.FASTQ)
        SnakemakeUtils.add_pickle_input(bwa_mem, 'INDEX_GENOME_PREFIX', input.FASTA_GENOME)
        SnakemakeUtils.add_pickle_input(sam_to_bam, 'FASTA_REF', input.FASTA_REF)

        bwa_mem.update_parameters(
            threads = threads,
            **config['rule_params']['alignment'][rule]
        )

        sam_to_bam.update_parameters(
            output_filename = params.output_file,
            threads = threads
        )

        pipeutils.run_as_pipe([bwa_mem, sam_to_bam], params.working_dir)

        SnakemakeUtils.dump_tool_output(sam_to_bam, "BAM", output.BAM)

rule picard_fastq_to_ubam:
    input:
        FASTQ = Path(config['working_dir']) / 'input' / "{input_basename}.fastq.gz.io",
    output:
        BAM_UNMAPPED = Path(config['working_dir']) / "alignment" / "fastq_to_ubam" / "{input_basename}.unmapped.bam.io"
    params:
        working_dir = Path(config['working_dir']) / "alignment" / "fastq_to_ubam",
        output_file = lambda wildcards: f'{wildcards.input_basename}.unmapped.bam'
    run:
        from camel.app.tools.picard.fastqtoubam import FastqTouBam

        Path(params.working_dir).mkdir(exist_ok=True)

        fastq_to_ubam = FastqTouBam(camel)

        SnakemakeUtils.add_pickle_input(fastq_to_ubam, 'FASTQ_PE', input.FASTQ)
        step = Step(rule, fastq_to_ubam, camel, params.working_dir, config)
        fastq_to_ubam.update_parameters(
            sample_name = config["sample"],
            readgroup_name = config["sample"] + ".rg",
            output = params.output_file,
            **config['rule_params']['alignment'][rule]
        )
        fastq_to_ubam.update_java_options('-mx20G -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dpicard.useLegacyParser=false')
        step.run_step()
        SnakemakeUtils.dump_tool_output(fastq_to_ubam, "BAM", output.BAM_UNMAPPED)

rule picard_merge_bam:
    input:
        BAM_UNMAPPED = rules.picard_fastq_to_ubam.output.BAM_UNMAPPED,
        SAM_ALIGNED = rules.bwa_aln_to_bam.output.BAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io"
    output:
        BAM = Path(config['working_dir']) / "alignment" / "merge_bam" / "{input_basename}merged.aligned.unsorted.bam.io",
    params:
        working_dir = Path(config['working_dir']) / 'alignment' / 'merge_bam',
        output_file = lambda wildcards: f'{wildcards.input_basename}§merged.aligned.unsorted.bam'
    run:
        from camel.app.tools.picard.mergebamalignment import MergeBamAlignment

        Path(params.working_dir).mkdir(exist_ok=True)

        merge_bam = MergeBamAlignment(camel)
        SnakemakeUtils.add_pickle_input(merge_bam, 'BAM_ALIGNED', input.SAM_ALIGNED)
        SnakemakeUtils.add_pickle_input(merge_bam, 'BAM_UNMAPPED', input.BAM_UNMAPPED)
        SnakemakeUtils.add_pickle_input(merge_bam, 'FASTA_REF', input.FASTA_REF)
        step = Step(rule, merge_bam, camel, params.working_dir, config)
        merge_bam.update_parameters(
            **config['rule_params']['alignment'][rule],
            output = params.output_file
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(merge_bam, "BAM", output.BAM)

rule picard_mark_duplicates:
    """
    Aggregate aligned+merged flowcell BAM files and mark duplicates
    We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
    to avoid having to spend time just merging BAM files.
    """
    input:
        BAM = expand(rules.bwa_aln_to_bam.output.BAM, input_basename = config['input_basenames'])
    output:
        BAM_DUPMARKED = Path(config['working_dir']) / 'alignment' / 'mark_duplicates' / f'{config["sample"]}.aligned.unsorted.duplicates_marked.bam.io',
        metrics = Path(config['working_dir']) / 'alignment' / 'mark_duplicates' / f"{config['sample']}.duplicate_metrics.txt.io"
    params:
        working_dir = Path(config['working_dir']) / 'alignment' / 'mark_duplicates',
        output_file = f'{config["sample"]}.aligned.unsorted.duplicates_marked.bam'
    run:
        from camel.app.tools.picard.markduplicates import MarkDuplicates

        mark_duplicates = MarkDuplicates(camel)
        mark_duplicates.add_input_files({"BAM": [SnakemakeUtils.load_object(path)[0] for path in input.BAM]})
        step = Step(rule, mark_duplicates, camel, params.working_dir, config)
        mark_duplicates.update_parameters(
            output = params.output_file,
            matrics_output = output.metrics,
            **config['rule_params']['alignment'][rule],
        )
        mark_duplicates.update_java_options("-mx100G -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dpicard.useLegacyParser=false")
        step.run_step()
        SnakemakeUtils.dump_tool_output(mark_duplicates, "BAM", output.BAM_DUPMARKED)
        SnakemakeUtils.dump_tool_output(mark_duplicates, "METRICS", output.metrics)

rule picard_sort_sam:
    input:
        BAM_DUPMARKED = rules.picard_mark_duplicates.output.BAM_DUPMARKED,
    output:
        BAM_SORTED = Path(config['working_dir']) / "alignment" / "sort_dupmarked" / "aligned.duplicate_marked.sorted.bam.io",
    params:
        working_dir = Path(config['working_dir']) / 'alignment' / 'sort_dupmarked'
    run:
        from camel.app.tools.picard.sortsam import SortSam

        sort_sam = SortSam(camel)
        SnakemakeUtils.add_pickle_input(sort_sam, 'BAM', input.BAM_DUPMARKED)
        step = Step(rule, sort_sam, camel, params.working_dir, config)
        sort_sam.update_parameters(
            output = "aligned.duplicate_marked.sorted.bam",
            **config['rule_params']['alignment'][rule],
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(sort_sam, "BAM", output.BAM_SORTED)

rule picard_add_readgroups:
    input:
        BAM = rules.picard_sort_sam.output.BAM_SORTED
    output:
        BAM = Path(config['working_dir']) / 'alignment' / 'add_readgroups' / 'aligned_dupmarked_sorted_rgadded.bam.io'
    params:
        working_dir = Path(config['working_dir']) / 'alignment' / 'add_readgroups'
    run:
        from camel.app.tools.picard.addorreplacereadgroups import AddOrReplaceReadGroups

        add_rg = AddOrReplaceReadGroups(camel)
        SnakemakeUtils.add_pickle_inputs(add_rg, input)
        step = Step(rule, add_rg, camel, params.working_dir, config)
        add_rg.update_parameters(
            output = 'aligned_dupmarked_sorted_rgadded.bam',
            RG_id = f'{config["sample"]}.rg',
            RG_sample_name = config['sample'],
            **config['rule_params']['alignment'][rule],
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(add_rg, "BAM", output.BAM)

rule picard_set_tags:
    input:
        BAM = rules.picard_add_readgroups.output.BAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io"
    output:
        BAM = Path(config['working_dir']) / 'alignment' / 'set_tags' / 'aligned_dupmarked_sorted_rgadded_settags.bam.io'
    params:
        working_dir = Path(config['working_dir']) / 'alignment' / 'set_tags',
    run:
        from camel.app.tools.picard.setnmmdanduqtags import SetNmMdAndUqTags

        set_tags = SetNmMdAndUqTags(camel)
        SnakemakeUtils.add_pickle_inputs(set_tags, input)
        step = Step(rule, set_tags, camel, params.working_dir, config)
        set_tags.update_parameters(
            output = "aligned_dupmarked_sorted_rgadded_settags.bam",
            **config['rule_params']['alignment'][rule]
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(set_tags, "BAM", output.BAM)

rule gatk4_baserecalibrator:
    input:
        BAM = rules.picard_set_tags.output.BAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
        VCF_KNOWN_SNPS = Path(config['working_dir']) / "ref_input" / "dbsnp_vcf.io",
        VCF_KNOWN_INDELS = Path(config['working_dir']) / "ref_input" / "known_indels_vcf.io",
        TXT_intervals = Path(config['intervals_location']) / "interval_{i}.intervals.io"
    output:
        TXT_RecalibrationTable = Path(config['working_dir']) / "alignment" / "bqsr" / "{i}_recal_data.csv.io"
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / "alignment" / "bqsr",
        output_file = lambda wildcards: f'{wildcards.i}_recal_data.csv'
    run:
        from camel.app.tools.gatk4.gatk4baserecalibrator import GATK4BaseRecalibrator

        bqsr = GATK4BaseRecalibrator(camel)
        SnakemakeUtils.add_pickle_inputs(bqsr, input)
        step = Step(rule, bqsr, camel, params.working_dir, config)
        bqsr.update_parameters(
            **config['rule_params']['alignment'][rule],
            output = params.output_file
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(bqsr, "TXT_RecalibrationTable", output.TXT_RecalibrationTable)

rule gatk4_gather_bqsr_reports:
    input:
        bqsr_report_interval = expand(rules.gatk4_baserecalibrator.output.TXT_RecalibrationTable, i = config["intervals"])
    output:
        bqsr_report_gathered = Path(config['working_dir']) / "alignment" / "gather_bqsr_reports" / "recal_data.csv.io",
    params:
        working_dir = Path(config['working_dir']) / "alignment" / "gather_bqsr_reports",
    run:
        from camel.app.tools.gatk4.gatk4gatherbqsrreports import GATK4GatherBQSRReports

        Path(params.working_dir).mkdir(exist_ok=True)

        gather_bqsr = GATK4GatherBQSRReports(camel)
        step = Step(rule, gather_bqsr, camel, params.working_dir, config)
        gather_bqsr.add_input_files({"TXT_intervals": [SnakemakeUtils.load_object(path)[0] for path in input.bqsr_report_interval]})
        gather_bqsr.update_parameters(
            output = "recal_data.csv"
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(gather_bqsr, 'TXT_RecalibrationTable', output.bqsr_report_gathered)

rule gatk4_apply_bqsr:
    input:
        BAM = rules.gatk4_baserecalibrator.input.BAM, # Ensure it gets same input as baserecalibrator
        BQSR = rules.gatk4_gather_bqsr_reports.output.bqsr_report_gathered,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
        TXT_intervals = Path(config['intervals_location']) / "interval_{i}.intervals.io"
    output:
        BAM = Path(config['working_dir']) / "alignment" / "apply_bqsr" / "{i}_sorted.bam.io"
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / "alignment" / "apply_bqsr",
        output_file = lambda wildcards: f'{wildcards.i}_sorted.bam'
    run:
        from camel.app.tools.gatk4.gatk4applybqsr import GATK4ApplyBQSR

        Path(params.working_dir).mkdir(exist_ok=True)

        apply_bqsr = GATK4ApplyBQSR(camel)
        SnakemakeUtils.add_pickle_inputs(apply_bqsr, input)
        step = Step(rule, apply_bqsr, camel, params.working_dir, config)
        apply_bqsr.update_parameters(
            **config['rule_params']['alignment'][rule],
            output = params.output_file
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(apply_bqsr, "BAM", output.BAM)

rule picard_gather_sorted_bam:
    input:
        bqsr_BAM_interval = expand(rules.gatk4_apply_bqsr.output.BAM, i = config["intervals"])
    output:
        BAM = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam" / "bqsr_gathered_sorted.bam.io",
    params:
        working_dir = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam"
    run:
        from camel.app.tools.picard.gatherbamfiles import GatherBamFiles

        Path(params.working_dir).mkdir(exist_ok=True)

        gather_bam = GatherBamFiles(camel)
        step = Step(rule, gather_bam, camel, params.working_dir, config)
        gather_bam.add_input_files({"BAMs": [SnakemakeUtils.load_object(path)[0] for path in input.bqsr_BAM_interval]})
        gather_bam.update_parameters(
            **config['rule_params']['alignment'][rule],
            output = "gathered_sorted.bam"
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(gather_bam, 'BAM', output.BAM)
