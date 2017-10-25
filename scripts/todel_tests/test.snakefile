import os
from yaml.error import Mark

from app.camel import Camel
from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue
from app.pipeline.snakestep import SnakeStep
from app.snakemake.snakemakeutils import SnakemakeUtils

camel = Camel()
working_dir = config['working_dir']

def prepare_addreadgroups_input(wildcards):
    """
    Prepares input for addreadgroups rule. Acts as a static fork in workflow based on config file (execution of markduplicates or not).
    :param wildcards: 
    :return: 
    """
    if config["run_markDuplicates"]:
        BAM = os.path.join(working_dir, "markduplicates/bam.io")
    elif not config["run_markDuplicates"]:
        BAM = os.path.join(working_dir, "picardsortbam/sortedbam.io")
    return BAM


rule all:
    # This rule makes sure that all other rules are executed.
    input:
        os.path.join(working_dir, "mutect1/txt.io")

rule prepare_initial_input:
    input:
        FASTQ=config['fastq_pe'],
        # DIR_OUT=[os.path.dirname(config['report'])]
    output:
        FASTQ=os.path.join(working_dir, "initial_input/fastq.io"),
        # DIR_OUT=os.path.join(working_dir, "initial_input/dir_html.io")
    run:
        SnakemakeUtils.pickle_snake_input(input, output)


rule prepare_bwa_input:
    input:
        FASTA_GENOME=config['ref_genome']
        # FASTA_GENOME=ToolIODb('broad_b37_human_Genome_1K_v37')
    output:
        FASTA_GENOME=os.path.join(working_dir, "initial_input/human_g1k_v37.fa.io")
    run:
        IO_FASTA_GENOME = [ToolIOValue(input.FASTA_GENOME)]
        SnakemakeUtils.dump_object(IO_FASTA_GENOME, output.FASTA_GENOME)


rule bwa_alignment:
    input:
        FASTQ=os.path.join(working_dir, "initial_input/fastq.io"),
        FASTA_GENOME=os.path.join(working_dir, "initial_input/human_g1k_v37.fa.io"),
    output:
        SAM=os.path.join(working_dir, "bwa_alignment/sam.io")
    threads: 16
    params:
        working_dir=os.path.join(working_dir, "bwa_alignment")
    run:
        from app.tools.bwa.bwamap import BWAMap
        bwa_mem = BWAMap(camel)
        SnakemakeUtils.add_pickle_input(bwa_mem, 'FASTQ_PE', input.FASTQ)
        SnakemakeUtils.add_pickle_input(bwa_mem, 'INDEX_GENOME_PREFIX', input.FASTA_GENOME)
        bwa_mem.update_parameters(threads=threads)
        # bwa_mem.run(os.path.join(working_dir, "bwa_alignment"))
        step = SnakeStep(rule, bwa_mem, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(bwa_mem, "SAM", output.SAM)


rule samtobam:
    input:
        SAM=os.path.join(working_dir, "bwa_alignment/sam.io"),
    output:
        BAM=os.path.join(working_dir, "bwa_alignment/bam.io"),
    run:
        from app.tools.samtools.samtoolsview import SamtoolsView
        smv = SamtoolsView(camel)
        SnakemakeUtils.add_pickle_input(smv,"SAM",input.SAM)
        smv.run(os.path.join(working_dir, "bwa_alignment"))
        SnakemakeUtils.dump_tool_output(smv, "BAM", output.BAM)


rule sortbam:
    input:
        BAM=os.path.join(working_dir, "bwa_alignment/bam.io"),
    output:
        BAM=os.path.join(working_dir, "bwa_alignment/sortedbam.io"),
    run:
        from app.tools.samtools.samtoolssort import SamtoolsSort
        sms = SamtoolsSort(camel)
        SnakemakeUtils.add_pickle_input(sms,"BAM",input.BAM)
        sms.run(os.path.join(working_dir, "bwa_alignment"))
        SnakemakeUtils.dump_tool_output(sms, "BAM", output.BAM)


rule indexbam:
    input:
        BAM=os.path.join(working_dir, "bwa_alignment/sortedbam.io"),
    output:
        BAM=os.path.join(working_dir, "samtools_index/bam_after_index.io"),
    run:
        from app.tools.samtools.samtoolsindex import SamtoolsIndex
        smi = SamtoolsIndex(camel)
        SnakemakeUtils.add_pickle_input(smi,"BAM",input.BAM)
        smi.run(os.path.join(working_dir, "samtools_index"))
        SnakemakeUtils.dump_tool_output(smi, "BAM", output.BAM)


rule picardsortbam:
    input:
        BAM=os.path.join(working_dir, "samtools_index/bam_after_index.io"),
    output:
        BAM=os.path.join(working_dir, "picardsortbam/sortedbam.io"),
    run:
        from app.tools.picard.sortsam import SortSam
        pss = SortSam(camel)
        pss.update_parameters(sort_order="coordinate")
        SnakemakeUtils.add_pickle_input(pss,"BAM",input.BAM)
        pss.run(os.path.join(working_dir, "picardsortbam"))
        SnakemakeUtils.dump_tool_output(pss, "BAM", output.BAM)


rule markduplicates:
    input:
        BAM=os.path.join(working_dir, "picardsortbam/sortedbam.io"),
    output:
        BAM=os.path.join(working_dir, "markduplicates/bam.io"),
        # METRICS=os.path.join(working_dir, "markduplicates/duplicates_mark.matrics.io"), not implemented in tool and maybe not useful. will see.
    run:
        from app.tools.picard.markduplicates import MarkDuplicates
        pmd = MarkDuplicates(camel)
        SnakemakeUtils.add_pickle_input(pmd, "BAM", input.BAM)
        pmd.run(os.path.join(working_dir, "markduplicates"))
        SnakemakeUtils.dump_tool_output(pmd, "BAM", output.BAM)
        # SnakemakeUtils.dump_tool_output(pss, "METRICS", output.METRICS)

rule addreadgroups:
    input:
        BAM=prepare_addreadgroups_input,
    output:
        BAM=os.path.join(working_dir, "addreadgroups/bam.io"),
    run:
        from app.tools.picard.addorreplacereadgroups import AddOrReplaceReadGroups
        parg = AddOrReplaceReadGroups(camel)
        parg.update_parameters(create_index='true')
        SnakemakeUtils.add_pickle_input(parg, "BAM", input.BAM)
        # parg.update_parameters(RG_sample_name=config["sample_name"])
        parg.run(os.path.join(working_dir, "addreadgroups"))
        SnakemakeUtils.dump_tool_output(parg, "BAM", output.BAM)

rule bamtobed:
    input:
        BAM=os.path.join(working_dir, "addreadgroups/bam.io"),
    output:
        BED=os.path.join(working_dir, "bamtobed/bed.io"),
    run:
        from app.tools.bedtools.bedtoolsbamtobed import BedtoolsBamToBed
        btb = BedtoolsBamToBed(camel)
        SnakemakeUtils.add_pickle_input(btb, "BAM", input.BAM)
        btb.run(os.path.join(working_dir, "bamtobed"))
        SnakemakeUtils.dump_tool_output(btb,"BED",output.BED)

rule generate_intervals:
    input:
        BED=os.path.join(working_dir, "bamtobed/bed.io"),
    output:
        BED=os.path.join(working_dir, "generate_intervals/bed.io"),
    run:
        from app.tools.bedtools.bedtoolsmerge import BedtoolsMerge
        btm = BedtoolsMerge(camel)
        SnakemakeUtils.add_pickle_input(btm, "BED", input.BED)
        btm.run(os.path.join(working_dir, "generate_intervals"))
        SnakemakeUtils.dump_tool_output(btm, "BED", output.BED)

rule realignertargetcreator:
    input:
        BAM=os.path.join(working_dir, "addreadgroups/bam.io"),
        BED=os.path.join(working_dir, "generate_intervals/bed.io"),
    output:
        INTERVALS=os.path.join(working_dir, "realignertargetcreator/intervals.io"),
    run:
        from app.tools.gatk.gatkrealignertargetcreator import GATKRealignerTargetCreator
        from app.io.tooliodb import ToolIODb
        grtc = GATKRealignerTargetCreator(camel)
        # add default human genome fasta file
        grtc.add_input_files({"FASTA_REF":[ToolIODb('broad_b37_human_Genome_1K_v37')]})
        SnakemakeUtils.add_pickle_input(grtc,"BAM",input.BAM)
        SnakemakeUtils.add_pickle_input(grtc,"TXT_intervals",input.BED)
        grtc.run(os.path.join(working_dir, "realignertargetcreator"))
        SnakemakeUtils.dump_tool_output(grtc, "TXT_realign_intervals", output.INTERVALS)

rule indelrealigner:
    input:
        INTERVALS=os.path.join(working_dir, "realignertargetcreator/intervals.io"),
        BAM=os.path.join(working_dir, "addreadgroups/bam.io"),
        BED=os.path.join(working_dir, "generate_intervals/bed.io"),
    output:
        BAM=os.path.join(working_dir, "indelrealigner/bam.io"),
    run:
        from app.tools.gatk.gatkindelrealigner import GATKIndelRealigner
        from app.io.tooliodb import ToolIODb
        gir=GATKIndelRealigner(camel)
        gir.add_input_files({"FASTA_REF":[ToolIODb('broad_b37_human_Genome_1K_v37')]})
        SnakemakeUtils.add_pickle_input(gir,"TXT_intervals",input.BED)
        SnakemakeUtils.add_pickle_input(gir,"BAM",input.BAM)
        SnakemakeUtils.add_pickle_input(gir,"TXT_realign_intervals",input.INTERVALS)
        gir.run(os.path.join(working_dir, "indelrealigner"))
        SnakemakeUtils.dump_tool_output(gir,"BAM",output.BAM)

rule basequalityrecalibration:
    input:
        BAM=os.path.join(working_dir, "indelrealigner/bam.io"),
        BED=os.path.join(working_dir, "generate_intervals/bed.io"),
    output:
        TXT=os.path.join(working_dir, "basequalityrecalibration/txt.io"),
    threads: 5
    run:
        from app.tools.gatk.gatkbaserecalibrator import GATKBaseRecalibrator
        from app.io.tooliodb import ToolIODb
        bqsr=GATKBaseRecalibrator(camel)
        bqsr.add_input_files({"FASTA_REF":[ToolIODb('broad_b37_human_Genome_1K_v37')]})
        bqsr.add_input_files({"VCF_KNOWN_SNPS":[ToolIODb('broad_b37_snps_high_confidence')]})
        bqsr.add_input_files({"VCF_KNOWN_INDELS":[ToolIODb('broad_b37_indels_gold_standard')]})
        bqsr.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(bqsr,"BAM",input.BAM)
        SnakemakeUtils.add_pickle_input(bqsr,"TXT_intervals", input.BED)
        bqsr.run(os.path.join(working_dir, "basequalityrecalibration"))
        SnakemakeUtils.dump_tool_output(bqsr,"TXT_RecalibrationTable",output.TXT)

rule printreads:
    input:
        BAM=os.path.join(working_dir, "indelrealigner/bam.io"),
        TXT=os.path.join(working_dir, "basequalityrecalibration/txt.io"),
    output:
        BAM=os.path.join(working_dir, "printreads/bam.io"),
    threads: 5
    run:
        from app.tools.gatk.gatkprintreads import GATKPrintReads
        from app.io.tooliodb import ToolIODb
        gpr=GATKPrintReads(camel)
        gpr.add_input_files({"FASTA_REF":[ToolIODb('broad_b37_human_Genome_1K_v37')]})
        gpr.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(gpr,"BAM",input.BAM)
        SnakemakeUtils.add_pickle_input(gpr,"BQSR",input.TXT)
        gpr.run(os.path.join(working_dir, "printreads"))
        SnakemakeUtils.dump_tool_output(gpr,"BAM",output.BAM)

rule basequalityrecalibration2:
    input:
        BAM = os.path.join(working_dir, "printreads/bam.io"),
        BED = os.path.join(working_dir, "generate_intervals/bed.io"),
    output:
        TXT = os.path.join(working_dir, "basequalityrecalibration2/txt.io"),
    threads: 5
    run:
        from app.tools.gatk.gatkbaserecalibrator import GATKBaseRecalibrator
        from app.io.tooliodb import ToolIODb

        bqsr = GATKBaseRecalibrator(camel)
        bqsr.add_input_files({"FASTA_REF": [ToolIODb('broad_b37_human_Genome_1K_v37')]})
        bqsr.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(bqsr, "BAM", input.BAM)
        SnakemakeUtils.add_pickle_input(bqsr, "TXT_intervals", input.BED)
        bqsr.run(os.path.join(working_dir, "basequalityrecalibration2"))
        SnakemakeUtils.dump_tool_output(bqsr, "TXT_RecalibrationTable", output.TXT)

rule analyzecovariates:
    input:
        TXT_BEFORE=os.path.join(working_dir, "basequalityrecalibration/txt.io"),
        TXT_AFTER=os.path.join(working_dir, "basequalityrecalibration2/txt.io"),
    output:
        PDF=os.path.join(working_dir, "analyzecovariates/pdf.io"),
    run:
        from app.tools.gatk.gatkanalyzecovariates import GATKAnalyzeCovariates

        gac = GATKAnalyzeCovariates(camel)
        SnakemakeUtils.add_pickle_input(gac,"TXT_TABLE_BEFORE",input.TXT_BEFORE)
        SnakemakeUtils.add_pickle_input(gac,"TXT_TABLE_AFTER",input.TXT_AFTER)
        gac.run(os.path.join(working_dir, "analyzecovariates"))
        SnakemakeUtils.dump_tool_output(gac, "PDF", output.PDF)

rule mutect1:
    input:
        BAM=os.path.join(working_dir, "printreads/bam.io"),
        BED = os.path.join(working_dir, "generate_intervals/bed.io"),
    output:
        TXT=os.path.join(working_dir, "mutect1/txt.io"),
        VCF=os.path.join(working_dir, "mutect1/vcf.io"),
    run:
        from app.tools.mutect.mutect1 import Mutect1
        from app.io.tooliodb import ToolIODb
        mut=Mutect1(camel)
        SnakemakeUtils.add_pickle_input(mut,'BAM_TUMOR',input.BAM)
        SnakemakeUtils.add_pickle_input(mut,"TXT_intervals",input.BED)
        mut.add_input_files({"FASTA_REF":[ToolIODb('broad_b37_human_Genome_1K_v37')]})
        mut.run(os.path.join(working_dir, "mutect1"))
        SnakemakeUtils.dump_tool_output(mut,'TXT_CALL_STATS',output.TXT)
        SnakemakeUtils.dump_tool_output(mut,'VCF',output.VCF)

