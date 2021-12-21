from pathlib import Path
#
from camel.resources.snakefile import trimming, trimming_illumina, trimming_iontorrent, assembly_spades, \
     quality_checks, contamination_check_kraken, gene_detection, pointfinder, variant_calling, variant_filtering, \
     sequence_typing
#from camel.scripts.salmonellapipeline.snakefile import serotype_detection

#######################
# Included Snakefiles #
#######################
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: trimming_iontorrent.SNAKEFILE_TRIMMING_IONTORRENT
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: quality_checks.SNAKEFILE_QUALITY_CHECKS

#########
# Rules #
#########
rule all:
    input:
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT

#include means that the smks are automatically executed if necessary; therefore i should now just ask for the output of the second quality check

rule select_fastq:
    """
    This rule creates an IO object with the trimmed FASTQ files.
    Other workflows such as Kraken or Assembly rely on this dictionary to get input files (PE or SE).
    """
    input:
        FASTQ_PE = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_DICT if trimming.get_read_type(config) == 'illumina' else [],
        FASTQ_SE = Path(config['working_dir']) / trimming_iontorrent.OUTPUT_TRIMMING_IONTORRENT_DICT if trimming.get_read_type(config) == 'iontorrent' else []
    output:
        IO_FASTQ = Path(config['working_dir']) / 'fq_dict.io'
    params:
        read_type = config['read_type']
    run:
        import shutil
        for key, fq in input.items():
            if len(fq) == 0:
                continue
            shutil.copyfile(fq, output.IO_FASTQ)
# rule all:
#     """
#     This rules ensures that the required output files are generated.
#     """
#
#     input:
#         'SRR17243381_1subset_fastqc.html',
#         'SRR17243381_2subset_fastqc.html'
#
#
# rule qualitycontrol1:
#     """
#     This rule will run FASTQC before trimming
#     """
#
#     input:
#         FQ_fwd = config['input']['fq_fwd'],
#         FQ_rev = config['input']['fq_rev']
#     output:
#         FORWARDFASTQC = 'SRR17243381_1subset_fastqc.html',
#         REVERSEFASTQC = 'SRR17243381_2subset_fastqc.html'
#     shell:
#         "fastqc {input.FQ_fwd} {input.FQ_rev};"

# rule trimming:
#     """
#     This rule will run Trimmomatic
#     """
#
#     input:
#
# rule qualitycontrol2:
#     """
#     This rule will run FASTQC after trimming
#     """
#
#     input:



#MK This rule should become active after trimming is implemented
# rule select_fastq:
#     """
#     This rule creates an IO object with the trimmed FASTQ files.
#     Other workflows such as Kraken or Assembly rely on this dictionary to get input files (PE or SE).
#     """
#     input:
#         FASTQ_PE = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_DICT if trimming.get_read_type(config) == 'illumina' else [],
#         FASTQ_SE = Path(config['working_dir']) / trimming_iontorrent.OUTPUT_TRIMMING_IONTORRENT_DICT if trimming.get_read_type(config) == 'iontorrent' else []
#     output:
#         IO_FASTQ = Path(config['working_dir']) / 'fq_dict.io'
#     params:
#         read_type = config['read_type']
#     run:
#         import shutil
#         for key, fq in input.items():
#             if len(fq) == 0:
#                 continue
#             shutil.copyfile(fq, output.IO_FASTQ)