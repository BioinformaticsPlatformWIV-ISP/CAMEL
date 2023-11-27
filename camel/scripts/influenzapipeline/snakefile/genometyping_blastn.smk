from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.tools.seqtk.seqtkconvert import SeqtkConvert
from camel.scripts.influenzapipeline.snakefile import assembly
from camel.scripts.influenzapipeline.snakefile import genometyping_blastn

camel = Camel.get_instance()



rule seqtk_subsample:
    """
    Runs the Seqtk subsampling on the data
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        FASTQ = Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_SUBSAMPLE_FASTQ,
        INFORMS = Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_SUBSAMPLE_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'genometyping' / 'seqtk_subsample',
        sample_name = config['sample_name']
    run:
        from camel.app.tools.seqtk.seqtksubsample import SeqtkSubsample
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.components.files.fastqutils import FastqUtils
        import logging

        fq_dict = SnakePipelineUtils.extracts_fq_input(Path(input.IO),key_pe='FASTQ_PE',key_se='FASTQ_SE')
        input_files = fq_dict['FASTQ_PE']

        base_count = 0
        for input_file in input_files:
            base_count += FastqUtils.count_bases(input_file.path)
        avg_depth = base_count / float(config['species_info']['genome_size'])
        target_depth = 200
        fraction = round(target_depth / avg_depth, 4)
        logger.info(f'Found a total of {base_count} bases in the input files. The average read depth based'
                     f'on this number is estimated to be {avg_depth}. Subsampling to ~200x (factor {fraction})')
        if fraction > 0.5:
            raise ValueError('Average read depth based on number of base pairs is less than 400 '
                             'which is not supported at this time!')

        subsample = SeqtkSubsample(camel)

        fq_dict = SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE', key_se='FASTQ_SE')
        input_files = fq_dict['FASTQ_PE']
        subsample.add_input_files({'FASTQ_PE': input_files})

        step = Step(str(rule), subsample, camel, params.running_dir)
        subsample.update_parameters(combine_output=True if config['analysis_type'] != 'assembly' else False,
            fraction=fraction)
        step.run_step()
        output_key = 'FASTQ' if config['analysis_type'] != 'assembly' else 'FASTQ_PE'
        SnakemakeUtils.dump_tool_output(subsample, output_key, Path(output.FASTQ))
        SnakemakeUtils.dump_tool_outputs(subsample, output, ['INFORMS'])

rule seqtk_convert:
    """
    Runs Seqtk convert on the data
    """
    input:
        FASTQ = Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_SUBSAMPLE_FASTQ
    output:
        FASTA = Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_CONVERT_FASTA,
        INFORMS = Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_CONVERT_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'genometyping' / 'seqtk_convert',
        sample_name = config['sample_name']
    run:
        def update_seqid(fasta):
            """
            Replaces space characters with underscores in the read id of the fasta
            :param fasta: fasta file whoes records to be updated
            :return: None
            """
            import shutil

            basename = Path(fasta).stem
            original_fasta = Path(fasta).parent / f'{basename}_original_ids.fa'
            shutil.move(fasta, original_fasta)
            with open(fasta, 'w') as outf:
                with open(original_fasta, 'r') as inf:
                    for l in inf:
                        if l.startswith('>'):
                            outf.write('_'.join(l.split(' ')))
                        else:
                            outf.write(l)

        convert = SeqtkConvert(camel)
        SnakemakeUtils.add_pickle_inputs(convert, input)
        step = Step(str(rule), convert, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(convert, output)
        update_seqid(convert.tool_outputs['FASTA'][0].path)

def get_blastn_query_input() -> Path:
    """
    Returns the input that is needed for the blastn genometyping rule (query sequences)
    :return: Path to the input
    """
    if config['analysis_type'] == 'alignment':
        return Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_CONVERT_FASTA
    elif config['analysis_type'] == 'assembly':
        io_path = Path(config['working_dir']) / 'genometyping_db.io'
        SnakemakeUtils.dump_object([ToolIOFile(Path(config['genometyping_db']))], io_path)
        return io_path
    else:
        raise ValueError('Invalid or no analysis type given in config file!')


def get_blastn_subject_input() -> Path:
    """
    Returns the input that is needed for the blastn genometyping rule (subject sequences)
    :return: Path to the input
    """
    if config['analysis_type'] == 'alignment':
        io_path = Path(config['working_dir']) / 'genometyping_db.io'
        SnakemakeUtils.dump_object([ToolIOFile(Path(config['genometyping_db']))], io_path)
        return io_path
    elif config['analysis_type'] == 'assembly':
        return Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA_BLAST_DB
    else:
        raise ValueError('Invalid or no analysis type given in config file!')

rule blastn_genometyping:
    """
    Runs Blastn on the data
    """
    input:
        FASTA = get_blastn_query_input(),
        DB_BLAST = get_blastn_subject_input()
    output:
        ASN = Path(config['working_dir']) / genometyping_blastn.OUTPUT_BLASTN_ASN,
        TSV = Path(config['working_dir']) / genometyping_blastn.OUTPUT_BLASTN_TSV,
        INFORMS = Path(config['working_dir']) / genometyping_blastn.OUTPUT_BLASTN_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'genometyping' / 'blastn'
    run:
        from camel.app.tools.blast.blastn import Blastn
        from camel.app.command.command import Command

        blastn = Blastn(camel)
        SnakemakeUtils.add_pickle_inputs(blastn, input)
        step = Step(str(rule), blastn, camel, params.running_dir)
        blastn.update_parameters(**config['rule_parameters']['blastn_genometyping'])
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blastn, output, keys=['ASN', 'INFORMS'])

        format_columns = ' '.join(['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
        output_file = params.running_dir / 'genometyping_blast.tsv'
        cmd = Command(f'module load blast; blast_formatter -archive {blastn.tool_outputs["ASN"][0].path} -outfmt "6 {format_columns}" > {output_file}')
        cmd.run(Path.cwd())
        SnakemakeUtils.dump_object([ToolIOFile(output_file)], Path(output.TSV))

rule blastn_genometyping_processing:
    """
    Processes the Blastn results
    """
    input:
        ASN = Path(config['working_dir']) / genometyping_blastn.OUTPUT_BLASTN_ASN
    output:
        INFORMS = Path(config['working_dir']) / genometyping_blastn.OUTPUT_BLASTN_PROCESSING_INFORMS,
        FASTA = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_FASTA_REF
    params:
        running_dir = Path(config['working_dir']) / 'genometyping' / 'blastn_processing'
    run:
        from camel.app.tools.pipelines.genome_typing.genometyping import GenomeTyping
        gt = GenomeTyping(camel)
        SnakemakeUtils.add_pickle_inputs(gt, input)
        ref_fasta = get_blastn_query_input() if config['analysis_type'] == 'assembly' else get_blastn_subject_input()
        SnakemakeUtils.add_pickle_inputs(gt, {'REF_FASTA': ref_fasta})
        # gt.add_input_files({'DB_BLAST': [ToolIOFile(config['genometyping_db'])]})
        step = Step(str(rule), gt, camel, params.running_dir)
        gt.update_parameters(**{'multi_segment': str(config['multi_segment']),
                                'seqIDParser_type': config['species_info']['seqIDParser_type'],
                                'genometyping_method': config['analysis_type'],
                                'genome_segments': config['species_info']['genome_segments'],
                                'random_seed': config['random_seed'],
                                'influenza_a': config['species_info']['name'] == 'Influenza A'})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(gt, output)


rule blastn_genometyping_reference_indexing:
    """
    Create indexes for the created reference genome (Bowtie2 and BWA)
    """
    input:
        FASTA_REF = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_FASTA_REF
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_INDEX_GENOME_PREFIX
    params:
        running_dir = Path(config['working_dir']) / 'genometyping' / 'blastn_processing'
    run:
        from camel.app.tools.bwa.bwaindex import BWAIndex
        from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex

        for indexer_class in [BWAIndex, Bowtie2Index]:
            indexer = indexer_class(camel)
            SnakemakeUtils.add_pickle_inputs(indexer, input)
            step = Step(str(rule), indexer, camel, params.running_dir)
            step.run_step()
            SnakemakeUtils.dump_tool_outputs(indexer, output)

        indexer = SamtoolsFastaIndex(camel)
        SnakemakeUtils.add_pickle_inputs(indexer, {'FASTA': input.FASTA_REF})
        step = Step(str(rule), indexer, camel, params.running_dir)
        step.run_step()

rule genometyping_report:
    """
    Creates the HTML report for the genome typing.
    """
    input:
        INFORMS_genometyping = Path(config['working_dir']) / genometyping_blastn.OUTPUT_BLASTN_PROCESSING_INFORMS,
        INFORMS_seqtksubsample = Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_SUBSAMPLE_INFORMS,
        TSV = Path(config['working_dir']) / genometyping_blastn.OUTPUT_BLASTN_TSV,
        FASTA = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_FASTA_REF
    output:
        VAL_HTML = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'genometyping' / 'report'
    run:
        from camel.app.tools.pipelines.genome_typing.reportergenometyping import ReporterGenomeTyping

        reporter = ReporterGenomeTyping(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule genometyping_export_summary:
    """
    Exports the summary information for the genome typing.
    """
    input:
        INFORMS_genometyping = Path(config['working_dir']) / genometyping_blastn.OUTPUT_BLASTN_PROCESSING_INFORMS,
        INFORMS_seqtksubsample = Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_SUBSAMPLE_INFORMS
    output:
        TSV = Path(config['working_dir'], genometyping_blastn.OUTPUT_GENOMETYPING_SUMMARY)
    run:
        gt_informs = SnakemakeUtils.load_object(Path(input.INFORMS_genometyping))
        seqtk_informs = SnakemakeUtils.load_object(Path(input.INFORMS_seqtksubsample))
        with open(output.TSV, 'w') as handle:
            handle.write(f'reads_for_subtyping\t{seqtk_informs["reads_count"]}\n')
            exp_segments = ','.join(gt_informs['expected_segments'])
            seg_recovered = ','.join(gt_informs['segment_coverage']['segment_covered'])
            seg_missing = ','.join(gt_informs['segment_coverage']['segment_missing'])
            seg_coverage = gt_informs['segment_coverage']['coverage']
            handle.write(f'segments_expected\t{exp_segments}\n')
            handle.write(f'segments_recovered\t{seg_recovered}\n')
            handle.write(f'segments_missing\t{seg_missing}\n')
            handle.write(f'segments_recovered_percentage\t{seg_coverage}\n')
            for segment in gt_informs['expected_segments']:
                if segment in gt_informs['segment_coverage']['segment_covered']:
                    handle.write(f'best_refseq_{segment}\t{gt_informs["segment_informs"][segment]["refseqid"]}\n')
                    candidates = ','.join(gt_informs['segment_informs'][segment]['candidates'])
                    handle.write(f'refseq_candidates_{segment}\t{candidates}\n')
            if 'hana_subtyping' in gt_informs:
                handle.write(f"subtype\t{gt_informs['hana_subtyping']['subtype']}\n")
                handle.write(f"subtype\t{gt_informs['hana_subtyping']['ha']}\n")
                handle.write(f"subtype\t{gt_informs['hana_subtyping']['na']}\n")