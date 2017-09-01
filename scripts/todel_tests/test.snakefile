import os
from yaml.error import Mark

from app.camel import Camel
from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue
from app.pipeline.snakestep import SnakeStep
from app.snakemake.snakemakeutils import SnakemakeUtils

camel = Camel()
working_dir = config['working_dir']

rule all:
    # This rule makes sure that all other rules are executed.
    input:
        os.path.join(working_dir, "addreadgroups/bam.io")


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
    run:
        from app.tools.bwa.bwamap import BWAMap
        bwa_mem = BWAMap(camel)
        SnakemakeUtils.add_pickle_input(bwa_mem, 'FASTQ_PE', input.FASTQ)
        SnakemakeUtils.add_pickle_input(bwa_mem, 'INDEX_GENOME_PREFIX', input.FASTA_GENOME)
        bwa_mem.update_parameters(threads=threads)
        bwa_mem.run(os.path.join(working_dir, "bwa_alignment"))
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

rule determine_addreadgroups_input:
    # choose input for addreadgroups.
    # Depends on execution of markDuplicates (see config file).
    # Implements a static fork in the workflow.
    output:
        ADDREADGROUPS_INPUT = os.path.join(working_dir, "addreadgroups/bam_input.io")
    run:
        if config["run_markDuplicates"]:
            BAM = [ToolIOFile(os.path.join(working_dir, "markduplicates/bam.io"))]
        elif not config["run_markDuplicates"]:
            BAM = [ToolIOFile(os.path.join(working_dir, "picardsortbam/sortedbam.io"))]
        SnakemakeUtils.dump_object(BAM, output.ADDREADGROUPS_INPUT)

rule addreadgroups:
    input:
        BAM=os.path.join(working_dir, "addreadgroups/bam_input.io")
    output:
        BAM=os.path.join(working_dir, "addreadgroups/bam.io")
    run:
        from app.tools.picard.addorreplacereadgroups import AddOrReplaceReadGroups
        parg = AddOrReplaceReadGroups(camel)
        SnakemakeUtils.add_pickle_input(parg, "BAM", input.BAM)
        parg.update_parameters(RG_sample_name=config["sample_name"])
        parg.run(os.path.join(working_dir, "addreadgroups"))
        SnakemakeUtils.dump_tool_output(parg, "BAM", output.BAM)


