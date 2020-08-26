from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.tools.seqtk.seqtkconvert import SeqtkConvert
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import trimming_illumina, deconseq
from camel.scripts.influenzapipeline.snakefile import subtyping_blastn
from camel.app.components.blasthit.influenzablastnasnparser import InfluenzaBlastnAsnParser


camel = Camel.get_instance()

rule seqtk_subsample:
    """
    Runs the Seqtk subsampling on the data
    """
    input:
        FASTQ_PE = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_CLEAN_PE if 'deconseq' in config['analyses'] else [],
        IO = Path(config['working_dir']) / 'fq_dict.io' if 'deconseq' not in config['analyses'] else []
    output:
        FASTQ = Path(config['working_dir']) / subtyping_blastn.OUTPUT_SEQTK_SUBSAMPLE_FASTQ,
        INFORMS = Path(config['working_dir']) / subtyping_blastn.OUTPUT_SEQTK_SUBSAMPLE_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'subtyping' / 'seqtk_subsample',
        sample_name = config['sample_name']
    run:
        from camel.app.tools.seqtk.seqtksubsample import SeqtkSubsample
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        try:
            fq_dict = SnakePipelineUtils.extracts_fq_input(input.IO, key_pe='FASTQ_PE', key_se='FASTQ_SE')
            input_file = fq_dict['FASTQ_PE'][0]
        except (AttributeError, TypeError):
            input_file = input.FASTQ_PE

        subsample = SeqtkSubsample(camel)
        subsample.add_input_files({'FASTQ': [input_file]})
        step = Step(rule, subsample, camel, params.running_dir, config)
        subsample.update_parameters(combine_output=True)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(subsample, output)

rule seqtk_convert:
    """
    Runs Seqtk convert on the data
    """
    input:
        FASTQ = Path(config['working_dir']) / subtyping_blastn.OUTPUT_SEQTK_SUBSAMPLE_FASTQ
    output:
        FASTA = Path(config['working_dir']) / subtyping_blastn.OUTPUT_SEQTK_CONVERT_FASTA,
        INFORMS = Path(config['working_dir']) / subtyping_blastn.OUTPUT_SEQTK_CONVERT_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'subtyping' / 'seqtk_convert',
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
        update_seqid(convert.get_outputs('FASTA')[0].path)

rule blastn_subtyping:
    """
    Runs Blastn on the data
    """
    input:
        FASTA = Path(config['working_dir']) / subtyping_blastn.OUTPUT_SEQTK_CONVERT_FASTA
    output:
        ASN = Path(config['working_dir']) / subtyping_blastn.OUTPUT_BLASTN_ASN,
        INFORMS = Path(config['working_dir']) / subtyping_blastn.OUTPUT_BLASTN_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'subtyping' / 'blastn'
    run:
        from camel.app.tools.blast.blastn import Blastn
        print(config)
        blastn = Blastn(camel)
        SnakemakeUtils.add_pickle_inputs(blastn, input)
        blastn.add_input_files({'DB_BLAST': [ToolIOFile(config['subtyping_db'])]})
        step = Step(rule, blastn, camel, params.running_dir, config)
        blastn.update_parameters(**config['rule_parameters']['blastn_subtyping'])
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(blastn, output)

rule blastn_subtyping_processing:
    """
    Processes the Blastn results
    """
    input:
        ASN = Path(config['working_dir']) / subtyping_blastn.OUTPUT_BLASTN_ASN
    output:
        INFORMS = Path(config['working_dir']) / subtyping_blastn.OUTPUT_BLASTN_PROCESSING_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'subtyping' / 'blastn_processing'
    run:
        from camel.app.tools.pipelines.segment_typing.segmenttyping import SegmentTyping
        st = SegmentTyping(camel)
        SnakemakeUtils.add_pickle_inputs(st, input)
        st.add_input_files({'DB_BLAST': [ToolIOFile(config['subtyping_db'])]})
        step = Step(rule, st, camel, params.running_dir, config)
        st.update_parameters(**{'multi_segment': str(config['multi_segment']),
                                'seqIDParser_type': config['species_info']['seqIDParser_type'],
                                'subtyping_method': 'blast',
                                'genome_segments': config['species_info']['genome_segments'],
                                'random_seed': config['random_seed']})
        import pickle
        pickle.dump(st, open('/scratch/rawinand/pickle_test', 'wb'))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(st, output)

        # inf = SnakemakeUtils.load_object(input.ASN)[0].path
        # inf_parser = InfluenzaBlastnAsnParser(inf, config['multi_segment'], config['species_info']['seqIDParser_type'], 'blast')
        # inf_parser.group_hits_per_segment()

