from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.tools.seqtk.seqtkconvert import SeqtkConvert
from camel.resources.snakefile import deconseq
from camel.scripts.influenzapipeline.snakefile import genometyping_blastn

camel = Camel.get_instance()

rule seqtk_subsample:
    """
    Runs the Seqtk subsampling on the data
    """
    input:
        FASTQ_PE = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_CLEAN_PE if 'deconseq' in config['analyses'] else [],
        IO = Path(config['working_dir']) / 'fq_dict.io' if 'deconseq' not in config['analyses'] else []
    output:
        FASTQ =Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_SUBSAMPLE_FASTQ,
        INFORMS =Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_SUBSAMPLE_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'genometyping' / 'seqtk_subsample',
        sample_name = config['sample_name']
    run:
        from camel.app.tools.seqtk.seqtksubsample import SeqtkSubsample
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        subsample = SeqtkSubsample(camel)

        try:
            fq_dict = SnakePipelineUtils.extracts_fq_input(input.IO, key_pe='FASTQ_PE', key_se='FASTQ_SE')
            input_files = fq_dict['FASTQ_PE']
            subsample.add_input_files({'FASTQ_PE': input_files})
        except (AttributeError, TypeError):
            SnakemakeUtils.add_pickle_inputs(subsample, input, keys=['FASTQ_PE'])

        step = Step(rule, subsample, camel, params.running_dir, config)
        subsample.update_parameters(combine_output=True)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(subsample, output)

rule seqtk_convert:
    """
    Runs Seqtk convert on the data
    """
    input:
        FASTQ =Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_SUBSAMPLE_FASTQ
    output:
        FASTA =Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_CONVERT_FASTA,
          INFORMS =Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_CONVERT_INFORMS
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
        step = Step(rule, convert, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(convert, output)
        update_seqid(convert.tool_outputs['FASTA'][0].path)

rule blastn_genometyping:
    """
    Runs Blastn on the data
    """
    input:
        FASTA = Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_CONVERT_FASTA
    output:
        ASN = Path(config['working_dir']) / genometyping_blastn.OUTPUT_BLASTN_ASN,
        TSV = Path(config['working_dir']) / genometyping_blastn.OUTPUT_BLASTN_TSV,
        INFORMS = Path(config['working_dir']) / genometyping_blastn.OUTPUT_BLASTN_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'genometyping' / 'blastn'
    run:
        from camel.app.tools.blast.blastn import Blastn
        from camel.app.command.command import Command
        import os

        blastn = Blastn(camel)
        SnakemakeUtils.add_pickle_inputs(blastn, input)
        blastn.add_input_files({'DB_BLAST': [ToolIOFile(config['genometyping_db'])]})
        step = Step(rule, blastn, camel, params.running_dir, config)
        blastn.update_parameters(**config['rule_parameters']['blastn_genometyping'])
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blastn, output, keys=['ASN', 'INFORMS'])

        format_columns = ' '.join(['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
        output_file = params.running_dir / 'genometyping_blast.tsv'
        cmd = Command(f'module load blast; blast_formatter -archive {blastn.tool_outputs["ASN"][0].path} -outfmt "6 {format_columns}" > {output_file}')
        cmd.run_command(os.getcwd())
        SnakemakeUtils.dump_object([ToolIOFile(output_file)], output.TSV)

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
        gt.add_input_files({'DB_BLAST': [ToolIOFile(config['genometyping_db'])]})
        step = Step(rule, gt, camel, params.running_dir, config)
        gt.update_parameters(**{'multi_segment': str(config['multi_segment']),
                                'seqIDParser_type': config['species_info']['seqIDParser_type'],
                                'genometyping_method': 'blast',
                                'genome_segments': config['species_info']['genome_segments'],
                                'random_seed': config['random_seed'],
                                'influenza_a': config['species_info']['name'] == 'Influenza A'})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(gt, output)


        # inf = SnakemakeUtils.load_object(input.ASN)[0].path
        # inf_parser = InfluenzaBlastnAsnParser(inf, config['multi_segment'], config['species_info']['seqIDParser_type'], 'blast')
        # inf_parser.group_hits_per_segment()

rule blastn_genometyping_reference_indexing:
    """
    Create indexes for the created reference genome (Bowtie2 and BWA)
    """
    input:
        FASTA_REF = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_FASTA_REF
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_INDEX_GENOME_PREFIX
    params:
        running_dir = Path(config['working_dir']) / 'genometyping' / 'ref_indexes'
    run:
        from camel.app.tools.bwa.bwaindex import BWAIndex
        from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index

        for indexer_class in [BWAIndex, Bowtie2Index]:
            indexer = indexer_class(camel)
            SnakemakeUtils.add_pickle_inputs(indexer, input)
            step = Step(rule, indexer, camel, params.running_dir, config)
            step.run_step()
            SnakemakeUtils.dump_tool_outputs(indexer, output)

rule genometyping_report:
    """
    Creates the HTML report for the genome typing.
    """
    input:
        INFORMS_genometyping = rules.blastn_genometyping_processing.output.INFORMS,
        INFORMS_seqtksubsample = rules.seqtk_subsample.output.INFORMS,
        TSV = rules.blastn_genometyping.output.TSV,
        FASTA = rules.blastn_genometyping_processing.output.FASTA
    output:
        VAL_HTML = Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'genometyping' / 'report'
    run:
        from camel.app.tools.pipelines.genome_typing.reportergenometyping import ReporterGenomeTyping

        reporter = ReporterGenomeTyping(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule genometyping_export_summary:
    """
    Exports the summary information for the genome typing.
    """
    input:
        INFORMS_genometyping = rules.blastn_genometyping_processing.output.INFORMS,
        INFORMS_seqtksubsample = rules.seqtk_subsample.output.INFORMS,
    output:
        TSV = Path(config['working_dir'], genometyping_blastn.OUTPUT_GENOMETYPING_SUMMARY)
    run:
        gt_informs = SnakemakeUtils.load_object(input.INFORMS_genometyping)
        seqtk_informs = SnakemakeUtils.load_object(input.INFORMS_seqtksubsample)
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