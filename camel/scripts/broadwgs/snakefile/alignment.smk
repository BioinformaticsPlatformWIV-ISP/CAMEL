from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

from camel.app.io.tooliofile import ToolIOFile

camel = Camel()

rule prepare_input:
    input:
        uBAM = Path(config['working_dir']) / "input" / config["sample"] / "{ubam}.unmapped.bam"
    output:
        uBAM = Path(config['working_dir']) / "input" / config["sample"]  / "{ubam}.unmapped.bam.io"
    run:
        io_uBAM = [ToolIOFile(input.uBAM)]
        SnakemakeUtils.dump_object(io_uBAM, str(output))

rule picard_samtofastq:
    input:
        BAM = Path(config['working_dir']) / "input" / config["sample"] / "{ubam}.unmapped.bam.io"
    output:
        FASTQ = Path(config['working_dir']) / "alignment" / "samtofastq" / "{ubam}.fastq.io",
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / "alignment" / "samtofastq" / wildcards.ubam
    run:
        from camel.app.tools.picard.samtofastq import SamToFastq

        Path(params.working_dir).mkdir(exist_ok=True)

        samtofastq = SamToFastq(camel)
        SnakemakeUtils.add_pickle_input(samtofastq,"BAM",input.BAM)
        step = Step(rule, samtofastq, camel, params.working_dir, config)
        samtofastq.update_parameters(interleave="true", non_pf="true")
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtofastq, output)

rule bwa_mapping:
    input:
        FASTQ = Path(config['working_dir']) / "alignment" / "samtofastq" / "{ubam}.fastq.io",
        FASTA_GENOME=Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value.io"
    output:
        SAM = Path(config['working_dir']) / "alignment" / "mapping" / "{ubam}.aligned.sam.io",
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / 'alignment' / "mapping" / wildcards.ubam,
        threads = 5,
    run:
        from camel.app.tools.bwa.bwamap import BWAMap

        Path(params.working_dir).mkdir(exist_ok=True)

        bwa_mem = BWAMap(camel)
        SnakemakeUtils.add_pickle_input(bwa_mem, 'FASTQ_INT', input.FASTQ)
        SnakemakeUtils.add_pickle_input(bwa_mem, 'INDEX_GENOME_PREFIX', input.FASTA_GENOME)
        step = Step(rule, bwa_mem, camel, params.working_dir, config)
        bwa_mem.update_parameters(threads=params.threads, verbose=3, interleaved='', deterministic_aln=100000000)
        step.run_step()
        SnakemakeUtils.dump_tool_output(bwa_mem, "SAM", output.SAM)

rule picard_merge_bam:
    input:
        SAM_ALIGNED = Path(config['working_dir']) / "alignment" / "mapping" / "{ubam}.aligned.sam.io",
        BAM_UNMAPPED = Path(config['working_dir']) / "input" / config["sample"] / "{ubam}.unmapped.bam.io",
        FASTA_GENOME_FILE = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io"
    output:
        BAM = Path(config['working_dir']) / "alignment" / "merge_bam" / "{ubam}.aligned.unsorted.io",
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / 'alignment' / 'merge_bam' / wildcards.ubam
    run:
        from camel.app.tools.picard.mergebamalignment import MergeBamAlignment

        Path(params.working_dir).mkdir(exist_ok=True)

        bam_aln = MergeBamAlignment(camel)
        SnakemakeUtils.add_pickle_input(bam_aln, 'BAM_ALIGNED', input.SAM_ALIGNED)
        SnakemakeUtils.add_pickle_input(bam_aln, 'BAM_UNMAPPED', input.BAM_UNMAPPED)
        SnakemakeUtils.add_pickle_input(bam_aln, 'FASTA_REF', input.FASTA_GENOME_FILE)
        step = Step(rule, bam_aln, camel, params.working_dir, config)
        bam_aln.update_parameters(
            validation_stringency = 'SILENT',
            expected_orientations = 'FR',
            attributes_to_retain = 'X0',
            attributes_to_remove = ' ATTRIBUTES_TO_REMOVE='.join(['NM', 'MD']), ## To do: update tool to take multiple inputs
            paired_run = 'true',
            sort_order = 'unsorted',
            is_bisulfite_sequence = 'false',
            aligned_reads_only = 'false',
            clip_adapters = 'false',
            max_records_in_ram = '2000000',
            add_mate_cigar = 'true',
            max_gaps = '-1',
            primary_alignment_strategy = 'MostDistant',
            unmapped_read_strategy = 'COPY_TO_TAG',
            aligner_paring = 'true',
            unmap_contaminant_reads = 'true',
            add_pg_tag_to_reads = 'false'
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(bam_aln, "BAM", output.BAM)


rule picard_mark_duplicates:
    input:
        BAM = expand(Path(config['working_dir']) / 'alignment' / 'merge_bam' / '{ubam}.aligned.unsorted.io', ubam = config['ubams'])
    output:
        BAM_DUPMARKED = Path(config['working_dir']) / 'alignment' / 'mark_duplicates' / (config["sample"] + 'aligned.unsorted.duplicates_marked.bam.io'),
        metrics = Path(config['working_dir']) / 'alignment' / 'metrics' / "duplicate_metrics.txt.io"
    params:
        working_dir = Path(config['working_dir']) / 'alignment' / 'mark_duplicates'
    run:
        from camel.app.tools.picard.markduplicates import MarkDuplicates

        mark_duplicates = MarkDuplicates(camel)
        mark_duplicates.add_input_files({"BAM": [SnakemakeUtils.load_object(path)[0] for path in input.BAM]})
        step = Step(rule, mark_duplicates, camel, params.working_dir, config)
        mark_duplicates.update_parameters(
            matrics_output = output.metrics,
            validation_stringency = 'SILENT',
            optical_duplicate_pixel_distance = '2500',
            assume_sort_order = "queryname",
            clear_dt = 'false',
            add_pg_tag_to_reads = 'false'
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(mark_duplicates, "BAM", output.BAM_DUPMARKED)
        SnakemakeUtils.dump_tool_output(mark_duplicates, "METRICS", output.metrics)

rule picard_sort_sam:
    input:
        BAM_DUPMARKED = Path(config['working_dir']) / 'alignment' / 'mark_duplicates' / (config["sample"] + 'aligned.unsorted.duplicates_marked.bam.io'),
    output:
        BAM_SORTED = Path(config['working_dir']) / "alignment" / "sort_dupmarked" / (config["sample"] + "aligned.duplicate_marked.sorted.bam.io"),
    params:
        working_dir = Path(config['working_dir']) / 'alignment' / 'sort_dupmarked'
    run:
        from camel.app.tools.picard.sortsam import SortSam

        sort_sam = SortSam(camel)
        SnakemakeUtils.add_pickle_input(sort_sam, 'BAM', input.BAM_DUPMARKED)
        step = Step(rule, sort_sam, camel, params.working_dir, config)
        sort_sam.update_parameters(
            sort_order = 'coordinate',
            create_index = 'true',
            create_md5_file = 'true',
            max_records_in_ram = 300000
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(sort_sam, "BAM", output.BAM_SORTED)

checkpoint create_sequence_grouping_intervals:
    input:
        DICT_GENOME = Path(config['working_dir']) / "ref_input" / "dictionary_genome_human.io",
    output:
        sequence_grouping = directory(Path(config['working_dir']) / "alignment"/ "interval_files" )
    run:
        Path(output.sequence_grouping).mkdir(exist_ok=True)

        dict_genome = SnakemakeUtils.load_object(input.DICT_GENOME)[0]

        with open(f"{dict_genome}", "r") as ref_dict_file:
            sequence_tuple_list = []
            longest_sequence = 0

            for line in ref_dict_file:
                if line.startswith("@SQ"):
                    line_split = line.split("\t")
                    # (Sequence_Name, Sequence_Length)
                    sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
            longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]

        # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
        # the last element after a :, so we add this as a sacrificial element.
        hg38_protection_tag = ":1+"

        # initialize the tsv string with the first sequence
        tsv_list = [(sequence_tuple_list[0][0] + hg38_protection_tag)]
        temp_size = sequence_tuple_list[0][1]
        tsv_print_list = []

        for sequence_tuple in sequence_tuple_list[1:]:
            if temp_size + sequence_tuple[1] <= longest_sequence:
                temp_size += sequence_tuple[1]
                tsv_list.append((sequence_tuple[0] + hg38_protection_tag))
            else:
                tsv_print_list.append(tsv_list)
                tsv_list = [(sequence_tuple[0] + hg38_protection_tag)]
                temp_size = sequence_tuple[1]
            tsv_list.append((sequence_tuple[0] + hg38_protection_tag))

        tsv_print_list.append(tsv_list)

        # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
        tsv_print_list.append(["unmapped"])

        for (n,interval) in enumerate(tsv_print_list):
            with open(Path(output.sequence_grouping) / ("interval_" + str(n) + ".intervals"), "w") as interval_file:
                interval_file.write("\n".join(interval))
                interval_file.close()

            interval_out = Path(output.sequence_grouping) / ("interval_" + str(n) + ".intervals")
            interval_io =  Path(output.sequence_grouping) / ("interval_" + str(n) + ".intervals.io")
            SnakemakeUtils.dump_object([ToolIOFile(interval_out)], str(interval_io))

rule gatk4_baserecalibrator:
    input:
        BAM = Path(config['working_dir']) / "alignment" / "sort_dupmarked" / (config["sample"] + "aligned.duplicate_marked.sorted.bam.io"),
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
        VCF_KNOWN_SNPS = Path(config['working_dir']) / "ref_input" / "dbsnp_vcf.io",
        TXT_intervals = Path(config['working_dir']) / "alignment" / "interval_files" / "interval_{i}.intervals.io"
    output:
        TXT_RecalibrationTable = Path(config['working_dir']) / "alignment" / "bqsr" / "{i}_recal.csv.io"
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / "alignment" / "bqsr" / wildcards.i
        #known_indels_sites_vcfs = references["known_indels_sites_vcfs"], ##TO DO
    run:
        from camel.app.tools.gatk4.gatk4baserecalibrator import GATK4BaseRecalibrator

        Path(params.working_dir).mkdir(exist_ok=True)
        bqsr = GATK4BaseRecalibrator(camel)
        SnakemakeUtils.add_pickle_inputs(bqsr, input)
        step = Step(rule, bqsr, camel, params.working_dir, config)
        bqsr.update_parameters(
            use_original_qualities = '',
            output = params.working_dir / 'recalibrationData.tabl'
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(bqsr, "TXT_RecalibrationTable", output.TXT_RecalibrationTable)

def get_intervals_bqsr(wildcards):
    dir_checkpoint = checkpoints.create_sequence_grouping_intervals.get(**wildcards).output.sequence_grouping
    all_intervals = glob_wildcards(os.path.join(dir_checkpoint, "interval_{i}.intervals.io")).i
    return expand(str(Path(config['working_dir']) / "alignment" / "bqsr" / "{i}_recal.csv.io"), i = all_intervals)

rule gatk4_gather_bqsr_reports:
    input:
        interval_files = get_intervals_bqsr
    output:
        recalibration_report_out = Path(config['working_dir']) / "alignment" / "gather_bqsr_reports" / "recal_data.csv.io",
    params:
        working_dir = Path(config['working_dir']) / "alignment" / "gather_bqsr_reports"
    run:
        # Workaround for snakemake issue #55
        interval_files_local = get_intervals_bqsr(wildcards)

        from camel.app.tools.gatk4.gatk4gatherbqsrreports import GATK4GatherBQSRReports

        Path(params.working_dir).mkdir(exist_ok=True)

        gather_bqsr = GATK4GatherBQSRReports(camel)
        step = Step(rule, gather_bqsr, camel, params.working_dir, config)
        gather_bqsr.add_input_files({"TXT_intervals": [SnakemakeUtils.load_object(path)[0] for path in interval_files_local]})
        step.run_step()
        SnakemakeUtils.dump_tool_output(gather_bqsr, 'TXT_RecalibrationTable', output.recalibration_report_out)

rule gatk4_apply_bqsr:
    input:
        BAM = Path(config['working_dir']) / "alignment" / "sort_dupmarked" / (config["sample"] + "aligned.duplicate_marked.sorted.bam.io"),
        BQSR = Path(config['working_dir']) / "alignment" / "gather_bqsr_reports" / "recal_data.csv.io",
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
        TXT_intervals = Path(config['working_dir']) / "alignment" / "interval_files" / "interval_{i}.intervals.io"
    output:
        BAM = Path(config['working_dir']) / "alignment" / "apply_bqsr" / (config["sample"] + "{i}_sorted_bam.io")
    params:
        working_dir = Path(config['working_dir']) / "alignment" / "apply_bqsr"
    run:
        from camel.app.tools.gatk4.gatk4applybqsr import GATK4ApplyBQSR

        Path(params.working_dir).mkdir(exist_ok=True)

        apply_bqsr = GATK4ApplyBQSR(camel)
        SnakemakeUtils.add_pickle_inputs(apply_bqsr, input)
        step = Step(rule, apply_bqsr, camel, params.working_dir, config)
        apply_bqsr.update_parameters(
            create_md5 = '',
            add_output_sam_PR = '',
            use_original_qualities = '',
            static_quantized_quals_multi = "10,20,30",
            output = params.working_dir / 'output.bam'
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(apply_bqsr, "BAM", output.BAM)

def get_intervals(wildcards):
    dir_checkpoint = checkpoints.create_sequence_grouping_intervals.get(**wildcards).output.sequence_grouping
    intervals = glob_wildcards(os.path.join(dir_checkpoint, "interval_{i}.intervals.io")).i
    return expand(str(Path(config['working_dir']) / "alignment" / "apply_bqsr" / (config["sample"] + "{i}_sorted_bam.io")), i = intervals)

rule picard_gather_sorted_bam:
    input:
        interval_files = get_intervals,
    output:
        BAM = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam" / (config["sample"] + ".bam.io"),
    params:
        working_dir = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam"
    run:
        # Workaround for snakemake issue #55
        interval_files_local = get_intervals(wildcards)

        from camel.app.tools.picard.gatherbamfiles import GatherBamFiles

        Path(params.working_dir).mkdir(exist_ok=True)

        gather_bam = GatherBamFiles(camel)
        step = Step(rule, gather_bam, camel, params.working_dir, config)
        gather_bam.add_input_files({"BAMs": [SnakemakeUtils.load_object(path)[0] for path in interval_files_local})
        gather_bam.update_parameters(
            create_index = "true",
            create_md5_file = "true"
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(gather_bam, 'BAM', output.BAM)
