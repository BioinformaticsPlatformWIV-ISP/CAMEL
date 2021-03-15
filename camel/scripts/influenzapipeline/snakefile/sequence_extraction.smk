from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import deconseq
from camel.scripts.influenzapipeline.snakefile import alignment
from camel.scripts.influenzapipeline.snakefile import genometyping_blastn
from camel.scripts.influenzapipeline.snakefile import sequence_extraction

camel = Camel.get_instance()

rule create_sequence_dictionary:
    input:
        FASTA_REF = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_FASTA_REF
    output:
        FASTA_REF = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_REF_SEQUENCE_DICTIONARY
    params:
        running_dir = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_INDEX_GENOME_PREFIX.parent
    threads: 6
    run:
        from camel.app.tools.picard.createsequencedictionary import CreateSequenceDictionary

        csd = CreateSequenceDictionary(camel)
        SnakemakeUtils.add_pickle_inputs(csd, input)
        step = Step(str(rule), csd, camel, params.running_dir,config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(csd, output)

rule add_or_replace_read_groups:
    """
    Runs Picard AddOrReplaceReadGroups
    """
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM if 'bam' not in config else config['bam']
    output:
        BAM = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_ADD_OR_REPLACE_BAM
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction' / 'add_or_replace'
    threads: 6
    run:
        from camel.app.tools.picard.addorreplacereadgroups import AddOrReplaceReadGroups

        add_or_replace = AddOrReplaceReadGroups(camel)
        SnakemakeUtils.add_pickle_inputs(add_or_replace, input)
        step = Step(str(rule), add_or_replace, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(add_or_replace, output)

rule fastq_to_sam:
    input:
        FASTQ_PE = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_CLEAN_PE if 'deconseq' in config.get('analyses', '') else config['FASTQ_PE'],
        IO = Path(config['working_dir']) / 'fq_dict.io' if 'deconseq' not in config.get('analyses', '') and 'FASTQ_PE' not in config else []
    output:
        BAM = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_FASTQ_TO_SAM
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction' / 'fastq_to_sam'
    threads: 6
    run:
        from camel.app.tools.picard.fastqtosam import FastqToSam

        fts = FastqToSam(camel)

        if input.IO:
            fq_dict = SnakePipelineUtils.extracts_fq_input(str(input.IO), key_pe='FASTQ_PE', key_se='FASTQ_SE')
            input_files = fq_dict['FASTQ_PE']
            fts.add_input_files({'FASTQ_PE': input_files})
        else:
            SnakemakeUtils.add_pickle_inputs(fts, input, keys=['FASTQ_PE'])

        step = Step(str(rule), fts, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fts, output)

rule merge_bam_alignment:
    input:
        FASTA_REF = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_REF_SEQUENCE_DICTIONARY if 'fasta_ref' not in config else config['fasta_ref'],
        BAM_ALIGNED = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_ADD_OR_REPLACE_BAM,
        BAM_UNMAPPED = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_FASTQ_TO_SAM
    output:
        BAM = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_MERGE_BAM_ALIGNMENT
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction' / 'merge_bam_alignment'
    threads: 6
    run:
        from camel.app.tools.picard.mergebamalignment import MergeBamAlignment

        mba = MergeBamAlignment(camel)
        SnakemakeUtils.add_pickle_inputs(mba, input)
        step = Step(str(rule), mba, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(mba, output)

rule mark_duplicates:
    input:
        BAM = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_MERGE_BAM_ALIGNMENT
    output:
        BAM = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_MARK_DUPLICATES
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction' / 'mark_duplicates'
    threads: 6
    run:
        from camel.app.tools.picard.markduplicates import MarkDuplicates

        md = MarkDuplicates(camel)
        SnakemakeUtils.add_pickle_inputs(md, input)
        step = Step(str(rule), md, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(md, output)

rule haplotypecaller:
    input:
        BAM = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_MARK_DUPLICATES,
        FASTA_REF = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_FASTA_REF if 'fasta_ref' not in config else config['fasta_ref']
    output:
        VCF = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_HAPLOTYPECALLER
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction' / 'haplotypecaller'
    threads: 6
    run:
        from camel.app.tools.gatk4.gatk4haplotypecaller import GATK4HaplotypeCaller

        hc = GATK4HaplotypeCaller(camel)
        SnakemakeUtils.add_pickle_inputs(hc, input)
        step = Step(str(rule), hc, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(hc, output)

rule variantfiltration:
    input:
        VCF = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_HAPLOTYPECALLER,
        FASTA_REF = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_FASTA_REF if 'fasta_ref' not in config else config['fasta_ref']
    output:
        VCF = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_VARIANTFILTRATION
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction' / 'variantfiltration'
    threads: 6
    run:
        from camel.app.tools.gatk4.gatk4variantfiltration import GATK4VariantFiltration

        vf = GATK4VariantFiltration(camel)
        SnakemakeUtils.add_pickle_inputs(vf, input)
        step = Step(str(rule), vf, camel, params.running_dir, config)
        vf.update_parameters(**{'filter-names': 'HardQC_filter,lowDP',
                                'filter-expressions': 'QD<2.0||FS>60.0||MQ<40.0||ReadPosRankSum<-8.0,DP<200'})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(vf, output)

rule selectvariants:
    input:
        VCF = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_VARIANTFILTRATION,
        FASTA_REF = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_FASTA_REF if 'fasta_ref' not in config else config['fasta_ref']
    output:
        VCF = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_SELECTVARIANTS
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction' / 'selectvariants'
    threads: 6
    run:
        from camel.app.tools.gatk4.gatk4selectvariants import GATK4SelectVariants

        sv = GATK4SelectVariants(camel)
        SnakemakeUtils.add_pickle_inputs(sv, input)
        step = Step(str(rule), sv, camel, params.running_dir, config)
        sv.update_parameters(**{'exclude-filtered': True})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sv, output)

rule se_samtools_depth:
    input:
        BAM = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_MARK_DUPLICATES
    output:
        TSV = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_SAMTOOLS_DEPTH
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction' / 'samtools_depth'
    threads: 6
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth

        sd = SamtoolsDepth(camel)
        SnakemakeUtils.add_pickle_inputs(sd, input)
        step = Step(str(rule), sd, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sd, output)

rule se_samtools_depth_stats:
    input:
        TXT = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_SAMTOOLS_DEPTH,
        FASTA_REF = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_FASTA_REF if 'fasta_ref' not in config else config['fasta_ref']
    output:
        INFORMS = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_SAMTOOLS_DEPTH_STATS
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction' / 'samtools_depth_stats'
    threads: 6
    run:
        from camel.app.tools.samtools.samtoolsdepthstatsanalyzer import SamtoolsDepthStatsAnalyzer

        sda = SamtoolsDepthStatsAnalyzer(camel)
        SnakemakeUtils.add_pickle_inputs(sda, input)
        step = Step(str(rule), sda, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sda, output)

rule vcf_indel_scan:
    input:
        VCF = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_SELECTVARIANTS
    output:
        INFORMS = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_VCF_INDEL_SCAN
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction' / 'vcf_indel_scan'
    threads: 6
    run:
        from camel.scripts.influenzapipeline.custom_tools.vcfindelscanner import VCFIndelScanner

        vcfis = VCFIndelScanner(camel)
        SnakemakeUtils.add_pickle_inputs(vcfis, input)
        step = Step(str(rule), vcfis, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(vcfis, output)

rule region_calculator:
    input:
        INFORMS_IndelScanner = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_VCF_INDEL_SCAN,
        INFORMS_DepthStats = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_SAMTOOLS_DEPTH_STATS
    output:
        TXT_intervals = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_REGION_CALCULATOR_INTERVALS,
        INFORMS = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_REGION_CALCULATOR_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction' / 'region_calculator'
    threads: 6
    run:
        from camel.scripts.influenzapipeline.custom_tools.alignmentextractionregioncalculator import AlignmentSeqExtractionRegionCalculator

        rc = AlignmentSeqExtractionRegionCalculator(camel)
        SnakemakeUtils.add_pickle_inputs(rc, input)
        step = Step(str(rule), rc, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(rc, output)

rule consensus_sequence:
    input:
        VCF = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_SELECTVARIANTS,
        TXT_Intervals = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_REGION_CALCULATOR_INTERVALS,
        FASTA_REF = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_FASTA_REF if 'fasta_ref' not in config else config['fasta_ref']
    output:
        FASTA = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction' / 'consensus_sequence'
    threads: 6
    run:
        from camel.app.tools.gatk4.gatk4fastaalternatereferencemaker import GATK4FastaAlternateReferenceMaker

        farm = GATK4FastaAlternateReferenceMaker(camel)
        SnakemakeUtils.add_pickle_inputs(farm, input)
        step = Step(str(rule), farm, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(farm, output)

rule consensus_sequence_index:
    """
    Create indexes for the created reference genome (Bowtie2 and BWA)
    """
    input:
        FASTA_REF = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE_INDEX_PREFIX
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction' / 'consensus_sequence'
    run:
        from camel.app.tools.bwa.bwaindex import BWAIndex
        from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index

        for indexer_class in [BWAIndex, Bowtie2Index]:
            indexer = indexer_class(camel)
            SnakemakeUtils.add_pickle_inputs(indexer, input)
            step = Step(str(rule), indexer, camel, params.running_dir, config)
            step.run_step()
            SnakemakeUtils.dump_tool_outputs(indexer, output)

rule iterative_consensus_generation:
    input:
        FASTA_REF = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE,
        BAM = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_FASTQ_TO_SAM,
        FASTQ_PE = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_CLEAN_PE if 'deconseq' in config.get('analyses', '') else config['FASTQ_PE'],
        IO = Path(config['working_dir']) / 'fq_dict.io' if 'deconseq' not in config.get('analyses', '') and 'FASTQ_PE' not in config else [],
        VCF = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_SELECTVARIANTS
    output:
        FASTA = Path(config['working_dir']) / sequence_extraction.OUTPUT_SEQ_EXTRACTION_CONSENSUS_SEQUENCE_ITERATIVE
    params:
        running_dir = Path(config['working_dir']) / 'seq_extraction_iterative'
    threads: 6
    run:
        import pickle
        import logging
        from camel.app.io.tooliofile import ToolIOFile
        from camel.scripts.influenzapipeline.consensuswrapper import ConsensusSequenceWrapper

        def reached_threshold(cw_current):
            if len(cw_current.output.vcf_diffs[0]) == 0 and len(cw_current.output.vcf_diffs[1]) == 0:
                logging.info('Finished iteration because there is no difference anymore between the last'
                             'two iterations')
                return True
            if run == 1:  # For the first run, no other runs to compare to are available
                return False
            # Check if the new run has the same variants as the run before (first check is for run 2
            # with cw_list = [run1, run2]).
            # Example with 3 variants more in runX vs runY and 1 variant more in runY vs runX
            # run1 vs run2: variants at pos 1, 5, 6 in run1 and at pos 8 in run2
            # run2 vs run3: variants at pos 8 in run 2 and at pos 1, 5, 6 in run3
            if cw_list[run-2].output.vcf_diffs[0] == cw_list[run-1].output.vcf_diffs[1] and \
                cw_list[run-2].output.vcf_diffs[1] == cw_list[run-1].output.vcf_diffs[0]:
                logging.info('Finished iteration because same variant differences were found'
                             'in the previous run')
                return True
            if len(cw_list) > 2:
                # Check if the same variant differences are returning from another iteration
                for i in range(1, len(cw_list)):
                    if cw_list[i-1].output.vcf_diffs[0] == cw_list[-1].output.vcf_diffs[0] and \
                        cw_list[i-1].output.vcf_diffs[1] == cw_list[-1].output.vcf_diffs[1]:
                        logging.info('Finished iteration because same variant differences were'
                                      'found as in a previous run')
                        return True
            if len(cw_list) > 15:
                raise ValueError('After more than 15 iterations, the consensus sequence has not converged!')
            return False

        if input.IO:
            fq_dict = SnakePipelineUtils.extracts_fq_input(str(input.IO), key_pe='FASTQ_PE', key_se='FASTQ_SE')
            input_fastq = fq_dict['FASTQ_PE']
        else:
            with open(input.FASTQ_PE, 'rb') as handle:
                input_fastq = pickle.load(handle)

        cw_list = []
        run = 1
        cw = ConsensusSequenceWrapper(Path(f'{params.running_dir}_{run}'))
        cw.run_workflow(input.FASTA_REF,input_fastq,input.VCF,config['aligner'])
        if not reached_threshold(cw):
            cw_list.append(cw)
            while not reached_threshold(cw):
                run += 1
                fasta_ref = cw.output.fasta_ref
                vcf = cw.output.vcf
                cw = ConsensusSequenceWrapper(Path(f'{params.running_dir}_{run}'))
                cw.run_workflow(fasta_ref,input_fastq,vcf,config['aligner'])
                cw_list.append(cw)
        SnakemakeUtils.dump_object([ToolIOFile(str(cw.output.fasta_ref))], output.FASTA)
