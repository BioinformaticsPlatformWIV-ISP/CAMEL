import os

from app.camel import Camel
from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue
from app.pipeline.snakestep import SnakeStep
from app.snakemake.snakemakeutils import SnakemakeUtils

camel = Camel()
working_dir = config['working_dir']


def prepare_addreadgroups_input(wildcards):
    """
    Prepares input for addreadgroups rule. 
    Acts as a static fork in workflow based on config file (execution of markduplicates or not). 
    :return: 
    """
    if config["run_markDuplicates"]:
        BAM = os.path.join(working_dir, "markduplicates/bam.io")
    else:
        BAM = os.path.join(working_dir, "picardsortbam/sortedbam.io")
    return BAM


rule all:
    """
    This rule makes sure that all other rules are executed.
    Last file to be generated in pipeline is vcf.io.
    """
    input:
        os.path.join(working_dir, "mutect1/vcf.io"),
        os.path.join(working_dir, "analyzecovariates/pdf.io")



rule prepare_initial_input:
    """ Prepare input for the pipeline: generates io file from PE or SE fastq file(s) """
    input:
        FASTQ=config['fastq'],
    output:
        FASTQ=os.path.join(working_dir, "initial_input/fastq.io"),
    run:
        SnakemakeUtils.pickle_snake_input(input, output)



rule prepare_references_io:
    """
    Prepare reference genome IO files for snakemake to use: generate io files for
    - the reference genome fasta (value and file objects)
    - the snp vcf file
    - the indel vcf file.
    Requires reference names from db_loc to be present in the config file.
    """
    output:
        FASTA_GENOME=os.path.join(working_dir, "initial_input/fasta_reference_human_value.io"),
        FASTA_GENOME_FILE=os.path.join(working_dir, "initial_input/fasta_reference_human.io"), # for tools that require ToolIOFile instead of ToolIOValue
        VCF_KNOWN_SNPS=os.path.join(working_dir, "initial_input/vcf_known_snps.io"),
        VCF_KNOWN_INDELS=os.path.join(working_dir, "initial_input/vcf_known_indels.io"),
    run:
        from app.io.tooliodb import ToolIODb
        FASTA_GENOME = config['fasta_ref'],
        VCF_KNOWN_SNPS = config['vcf_known_snps'],
        VCF_KNOWN_INDELS = config['vcf_known_indels'],
        IO_FASTA_GENOME = [ToolIOValue(ToolIODb(FASTA_GENOME).path)]
        IO_FILE_FASTA_GENOME = [ToolIOFile(ToolIODb(FASTA_GENOME).path)]
        IO_VCF_KNOWN_SNPS = [ToolIOFile(ToolIODb(VCF_KNOWN_SNPS).path)]
        IO_VCF_KNOWN_INDELS = [ToolIOFile(ToolIODb(VCF_KNOWN_INDELS).path)]
        SnakemakeUtils.dump_object(IO_FASTA_GENOME, output.FASTA_GENOME)
        SnakemakeUtils.dump_object(IO_FILE_FASTA_GENOME, output.FASTA_GENOME_FILE)
        SnakemakeUtils.dump_object(IO_VCF_KNOWN_SNPS, output.VCF_KNOWN_SNPS)
        SnakemakeUtils.dump_object(IO_VCF_KNOWN_INDELS, output.VCF_KNOWN_INDELS)


rule bwa_alignment:
    """
    Reads alignment using bwa mem.
    """
    input:
        FASTQ=os.path.join(working_dir, "initial_input/fastq.io"),
        FASTA_GENOME=os.path.join(working_dir, "initial_input/fasta_reference_human_value.io"),
    output:
        SAM=os.path.join(working_dir, "bwa_alignment/sam.io")
    threads: 16
    params:
        working_dir=os.path.join(working_dir, "bwa_alignment"),
    run:
        from app.tools.bwa.bwamap import BWAMap
        bwa_mem = BWAMap(camel)
        if config['PE']:
            SnakemakeUtils.add_pickle_input(bwa_mem, 'FASTQ_PE', input.FASTQ)
        if config['SE']:
            SnakemakeUtils.add_pickle_input(bwa_mem, 'FASTQ_SE', input.FASTQ)
        SnakemakeUtils.add_pickle_input(bwa_mem, 'INDEX_GENOME_PREFIX', input.FASTA_GENOME)
        bwa_mem.update_parameters(threads=threads)
        step = SnakeStep(rule, bwa_mem, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(bwa_mem, "SAM", output.SAM)


rule samtobam:
    """
    Sam-to-bam conversion (samtools).
    """
    input:
        SAM=os.path.join(working_dir, "bwa_alignment/sam.io"),
    output:
        BAM=os.path.join(working_dir, "bwa_alignment/bam.io"),
    params:
        working_dir = os.path.join(working_dir, "bwa_alignment")
    run:
        from app.tools.samtools.samtoolsview import SamtoolsView
        smv = SamtoolsView(camel)
        SnakemakeUtils.add_pickle_input(smv,"SAM",input.SAM)
        step = SnakeStep(rule, smv, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(smv, "BAM", output.BAM)


rule sortbam:
    """
    pre-indexing Bam sorting (samtools). 
    """
    input:
        BAM=os.path.join(working_dir, "bwa_alignment/bam.io"),
    output:
        BAM=os.path.join(working_dir, "bwa_alignment/sortedbam.io"),
    params:
        working_dir = os.path.join(working_dir, "bwa_alignment")
    run:
        from app.tools.samtools.samtoolssort import SamtoolsSort
        sms = SamtoolsSort(camel)
        if 'bam_output' in config:
            sms.update_parameters(output_filename = config['bam_output'])
        SnakemakeUtils.add_pickle_input(sms,"BAM",input.BAM)
        step = SnakeStep(rule, sms, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(sms, "BAM", output.BAM)


rule indexbam:
    """
    Bam indexing (samtools).
    """
    input:
        BAM=os.path.join(working_dir, "bwa_alignment/sortedbam.io"),
    output:
        BAM=os.path.join(working_dir, "samtools_index/bam_after_index.io"),
    params:
        working_dir = os.path.join(working_dir, "samtools_index")
    run:
        from app.tools.samtools.samtoolsindex import SamtoolsIndex
        smi = SamtoolsIndex(camel)
        SnakemakeUtils.add_pickle_input(smi,"BAM",input.BAM)
        step = SnakeStep(rule, smi, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(smi, "BAM", output.BAM)


rule picardsortbam:
    """
    Bam sorting (picard).
    """
    input:
        BAM=os.path.join(working_dir, "samtools_index/bam_after_index.io"),
    output:
        BAM=os.path.join(working_dir, "picardsortbam/sortedbam.io"),
    params:
        working_dir = os.path.join(working_dir, "picardsortbam")
    run:
        from app.tools.picard.sortsam import SortSam
        pss = SortSam(camel)
        SnakemakeUtils.add_pickle_input(pss,"BAM",input.BAM)
        step = SnakeStep(rule, pss, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(pss, "BAM", output.BAM)


rule markduplicates:
    """
    Optional reads deduplication (picard).
    """
    input:
        BAM=os.path.join(working_dir, "picardsortbam/sortedbam.io"),
    output:
        BAM=os.path.join(working_dir, "markduplicates/bam.io"),
    params:
        working_dir = os.path.join(working_dir, "markduplicates")
    run:
        from app.tools.picard.markduplicates import MarkDuplicates
        pmd = MarkDuplicates(camel)
        SnakemakeUtils.add_pickle_input(pmd, "BAM", input.BAM)
        step = SnakeStep(rule, pmd, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(pmd, "BAM", output.BAM)


rule addreadgroups:
    """
    Read group adding (picard)
    """
    input:
        BAM=prepare_addreadgroups_input,
    output:
        BAM=os.path.join(working_dir, "addreadgroups/bam.io"),
    params:
        working_dir = os.path.join(working_dir, "addreadgroups")
    run:
        from app.tools.picard.addorreplacereadgroups import AddOrReplaceReadGroups
        parg = AddOrReplaceReadGroups(camel)
        SnakemakeUtils.add_pickle_input(parg, "BAM", input.BAM)
        step = SnakeStep(rule, parg, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(parg, "BAM", output.BAM)


rule bamtobed:
    """
    Bam-to-bed conversion for interval generation for optimisation (bedtools).
    """
    input:
        BAM=os.path.join(working_dir, "addreadgroups/bam.io"),
    output:
        BED=os.path.join(working_dir, "bamtobed/bed.io"),
    params:
        working_dir = os.path.join(working_dir, "bamtobed")
    run:
        from app.tools.bedtools.bedtoolsbamtobed import BedtoolsBamToBed
        btb = BedtoolsBamToBed(camel)
        SnakemakeUtils.add_pickle_input(btb, "BAM", input.BAM)
        step = SnakeStep(rule, btb, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(btb,"BED",output.BED)


rule generate_intervals:
    """
    Interval generation for optimisation (bedtools).
    """
    input:
        BED=os.path.join(working_dir, "bamtobed/bed.io"),
    output:
        BED=os.path.join(working_dir, "generate_intervals/bed.io"),
    params:
        working_dir = os.path.join(working_dir, "generate_intervals")
    run:
        from app.tools.bedtools.bedtoolsmerge import BedtoolsMerge
        btm = BedtoolsMerge(camel)
        SnakemakeUtils.add_pickle_input(btm, "BED", input.BED)
        step = SnakeStep(rule, btm, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(btm, "BED", output.BED)


rule realignertargetcreator:
    """
    Indel realigner intervals creation (GATK).
    """
    input:
        BAM=os.path.join(working_dir, "addreadgroups/bam.io"),
        BED=os.path.join(working_dir, "generate_intervals/bed.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),

    output:
        INTERVALS=os.path.join(working_dir, "realignertargetcreator/intervals.io"),
    params:
        working_dir = os.path.join(working_dir, "realignertargetcreator")
    run:
        from app.tools.gatk.gatkrealignertargetcreator import GATKRealignerTargetCreator
        grtc = GATKRealignerTargetCreator(camel)
        SnakemakeUtils.add_pickle_input(grtc,"FASTA_REF",input.FASTA_REF)
        SnakemakeUtils.add_pickle_input(grtc,"BAM",input.BAM)
        SnakemakeUtils.add_pickle_input(grtc,"TXT_intervals",input.BED)
        step = SnakeStep(rule, grtc, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(grtc, "TXT_realign_intervals", output.INTERVALS)


rule indelrealigner:
    """
    Indel realignment (GATK).
    """
    input:
        INTERVALS=os.path.join(working_dir, "realignertargetcreator/intervals.io"),
        BAM=os.path.join(working_dir, "addreadgroups/bam.io"),
        BED=os.path.join(working_dir, "generate_intervals/bed.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),
    output:
        BAM=os.path.join(working_dir, "indelrealigner/bam.io"),
    params:
        working_dir = os.path.join(working_dir, "indelrealigner"),
    run:
        from app.tools.gatk.gatkindelrealigner import GATKIndelRealigner
        gir=GATKIndelRealigner(camel)
        SnakemakeUtils.add_pickle_input(gir,"FASTA_REF",input.FASTA_REF)
        SnakemakeUtils.add_pickle_input(gir,"TXT_intervals",input.BED)
        SnakemakeUtils.add_pickle_input(gir,"BAM",input.BAM)
        SnakemakeUtils.add_pickle_input(gir,"TXT_realign_intervals",input.INTERVALS)
        step = SnakeStep(rule, gir, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(gir,"BAM",output.BAM)


rule basequalityrecalibration:
    """
    Base quality recalibration (GATK).
    """
    input:
        BAM=os.path.join(working_dir, "indelrealigner/bam.io"),
        BED=os.path.join(working_dir, "generate_intervals/bed.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),
        VCF_KNOWN_SNPS=os.path.join(working_dir, "initial_input/vcf_known_snps.io"),
        VCF_KNOWN_INDELS=os.path.join(working_dir, "initial_input/vcf_known_indels.io"),
    output:
        TXT=os.path.join(working_dir, "basequalityrecalibration/txt.io"),
    threads: 5
    params:
        working_dir = os.path.join(working_dir, "basequalityrecalibration"),
    run:
        from app.tools.gatk.gatkbaserecalibrator import GATKBaseRecalibrator
        bqsr=GATKBaseRecalibrator(camel)
        SnakemakeUtils.add_pickle_input(bqsr,"FASTA_REF",input.FASTA_REF)
        SnakemakeUtils.add_pickle_input(bqsr,"VCF_KNOWN_SNPS",input.VCF_KNOWN_SNPS)
        SnakemakeUtils.add_pickle_input(bqsr,"VCF_KNOWN_INDELS",input.VCF_KNOWN_INDELS)
        SnakemakeUtils.add_pickle_input(bqsr,"BAM",input.BAM)
        SnakemakeUtils.add_pickle_input(bqsr,"TXT_intervals", input.BED)
        bqsr.update_parameters(threads=threads)
        step = SnakeStep(rule, bqsr, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(bqsr,"TXT_RecalibrationTable",output.TXT)


rule printreads:
    """
    Post-recalibration reads printing (GATK). 
    """
    input:
        BAM=os.path.join(working_dir, "indelrealigner/bam.io"),
        TXT=os.path.join(working_dir, "basequalityrecalibration/txt.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),
    output:
        BAM=os.path.join(working_dir, "printreads/bam.io"),
    threads: 5
    params:
        working_dir = os.path.join(working_dir, "printreads"),
    run:
        from app.tools.gatk.gatkprintreads import GATKPrintReads
        gpr=GATKPrintReads(camel)
        gpr.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(gpr,"FASTA_REF",input.FASTA_REF)
        SnakemakeUtils.add_pickle_input(gpr,"BAM",input.BAM)
        SnakemakeUtils.add_pickle_input(gpr,"BQSR",input.TXT)
        step = SnakeStep(rule, gpr, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(gpr,"BAM",output.BAM)


rule basequalityrecalibration2:
    """
    Second base quality recalibration for quality control (GATK).
    """
    input:
        BAM = os.path.join(working_dir, "printreads/bam.io"),
        BED = os.path.join(working_dir, "generate_intervals/bed.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),
        VCF_KNOWN_SNPS=os.path.join(working_dir, "initial_input/vcf_known_snps.io"),
        VCF_KNOWN_INDELS=os.path.join(working_dir, "initial_input/vcf_known_indels.io"),
    output:
        TXT = os.path.join(working_dir, "basequalityrecalibration2/txt.io"),
    threads: 5
    params:
        working_dir = os.path.join(working_dir, "basequalityrecalibration2"),
    run:
        from app.tools.gatk.gatkbaserecalibrator import GATKBaseRecalibrator
        bqsr = GATKBaseRecalibrator(camel)
        bqsr.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(bqsr,"FASTA_REF",input.FASTA_REF)
        SnakemakeUtils.add_pickle_input(bqsr,"VCF_KNOWN_SNPS",input.VCF_KNOWN_SNPS)
        SnakemakeUtils.add_pickle_input(bqsr,"VCF_KNOWN_INDELS",input.VCF_KNOWN_INDELS)
        SnakemakeUtils.add_pickle_input(bqsr,"BAM",input.BAM)
        SnakemakeUtils.add_pickle_input(bqsr,"TXT_intervals", input.BED)
        step = SnakeStep(rule, bqsr, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(bqsr, "TXT_RecalibrationTable", output.TXT)


rule analyzecovariates:
    """
    Covariates analysis and BQSR report generation (GATK) 
    """
    input:
        TXT_BEFORE=os.path.join(working_dir, "basequalityrecalibration/txt.io"),
        TXT_AFTER=os.path.join(working_dir, "basequalityrecalibration2/txt.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),
    output:
        PDF=os.path.join(working_dir, "analyzecovariates/pdf.io"),
    params:
        working_dir = os.path.join(working_dir, "analyzecovariates")
    run:
        from app.tools.gatk.gatkanalyzecovariates import GATKAnalyzeCovariates
        gac = GATKAnalyzeCovariates(camel)
        if 'covar_output' in config:
            gac.update_parameters(pdf_output = config['covar_output'])
        SnakemakeUtils.add_pickle_input(gac,"FASTA_REF",input.FASTA_REF)
        SnakemakeUtils.add_pickle_input(gac,"TXT_TABLE_BEFORE",input.TXT_BEFORE)
        SnakemakeUtils.add_pickle_input(gac,"TXT_TABLE_AFTER",input.TXT_AFTER)
        step = SnakeStep(rule, gac, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(gac, "PDF", output.PDF)


rule mutect1:
    """
    Variant calling (mutect1).
    """
    input:
        BAM=os.path.join(working_dir, "printreads/bam.io"),
        BED = os.path.join(working_dir, "generate_intervals/bed.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),
        VCF_KNOWN_SNPS=os.path.join(working_dir, "initial_input/vcf_known_snps.io"),
    output:
        TXT=os.path.join(working_dir, "mutect1/txt.io"),
        VCF=os.path.join(working_dir, "mutect1/vcf.io"),
    params:
        working_dir = os.path.join(working_dir, "mutect1"),
    run:
        from app.tools.mutect.mutect1 import Mutect1
        mut=Mutect1(camel)
        SnakemakeUtils.add_pickle_input(mut,'BAM_TUMOR',input.BAM)
        SnakemakeUtils.add_pickle_input(mut,"TXT_intervals",input.BED)
        SnakemakeUtils.add_pickle_input(mut, "FASTA_REF", input.FASTA_REF)
        SnakemakeUtils.add_pickle_input(mut, "VCF_DBSNP", input.VCF_KNOWN_SNPS)
        if 'txt_output' in config:
            mut.update_parameters(output_callstats_file = config['txt_output'])
        if 'vcf_output' in config:
            mut.update_parameters(output_vcf_file = config['vcf_output'])
        if 'downsampling_target' in config:
            mut.update_parameters(downsampling_coverage_target=config['downsampling_target'])
        if 'downsampling_type' in config:
            mut.update_parameters(downsampling_type=config['downsampling_type'])
        step = SnakeStep(rule, mut, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(mut,'TXT_CALL_STATS',output.TXT)
        SnakemakeUtils.dump_tool_output(mut,'VCF',output.VCF)
