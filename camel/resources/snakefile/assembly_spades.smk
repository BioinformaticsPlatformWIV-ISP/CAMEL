from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import assembly_spades

camel = Camel.get_instance()


rule assembly_spades_run:
    """
    De-novo assembly using SPAdes.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        FASTA_Contig = Path(config['working_dir']) / 'assembly_spades' / 'spades' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_INFORMS
    params:
        dir_ = Path(config['working_dir']) / 'assembly_spades' / 'spades',
        spades_options = config.get('assembly', {}).get('spades', {}),
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE'
    threads: 8
    priority: 1
    run:
        from camel.app.tools.spades.spades import SPAdes
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        spades = SPAdes(camel)

        # Reformat FASTQ dictionary
        fq_dict = SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE_1', keys_se=[
            'FASTQ_SE_1', 'FASTQ_SE_2'], key_se='FASTQ_SE_1', drop_empty=True, read_type=params.read_type)
        spades.add_input_files(fq_dict)
        step = Step(str(rule), spades, camel, params.dir_)
        spades.update_parameters(**params.spades_options)
        spades.update_parameters(isolate=True, careful=False, threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spades, output)

rule assembly_spades_filter_contig_length:
    """
    Filters out the small contigs.
    """
    input:
        FASTA = rules.assembly_spades_run.output.FASTA_Contig
    output:
        FASTA = Path(config['working_dir']) / 'assembly_spades' / 'filtering' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FILTERING_INFORMS
    params:
        dir_ = Path(config['working_dir']) / 'assembly_spades' / 'filtering',
        min_contig_length = config.get('assembly', {}).get('min_contig_length', 0)
    run:
        from camel.app.tools.seqtk.seqtkseq import SeqtkSeq
        seqtk = SeqtkSeq(camel)
        SnakemakeUtils.add_pickle_inputs(seqtk, input)
        step = Step(str(rule), seqtk, camel, params.dir_)
        seqtk.update_parameters(output_filename='assembly_filtered.fasta', min_length=params.min_contig_length)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqtk, output)

rule assembly_spades_bt2_index:
    """
    Creates a bowtie2 index for the assembly.
    """
    input:
        FASTA_REF = rules.assembly_spades_filter_contig_length.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'assembly_spades' / 'bowtie2' / 'genome_prefix.io'
    params:
        dir_ = Path(config['working_dir']) / 'assembly_spades' / 'bowtie2'
    run:
        from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index
        bowtie2_index = Bowtie2Index(camel)
        step = Step(str(rule), bowtie2_index, camel, params.dir_)
        SnakemakeUtils.add_pickle_inputs(bowtie2_index, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bowtie2_index, output)

rule assembly_spades_bt2_map:
    """
    Maps the reads against the assembled contigs.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
        INDEX_GENOME_PREFIX = rules.assembly_spades_bt2_index.output.INDEX_GENOME_PREFIX
    output:
        BAM = Path(config['working_dir']) / 'assembly_spades' / 'bowtie2' / 'bam.io',
        INFORMS = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_MAPPING_INFORMS
    params:
        dir_ = Path(config['working_dir']) / 'assembly_spades' / 'bowtie2',
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE'
    threads: 8
    run:
        from camel.app.components.pipelines import pipeutils
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        # Bowtie 2
        bowtie2_map = Bowtie2Map(camel)
        bowtie2_map.add_input_files(SnakePipelineUtils.extracts_fq_input(
            Path(input.IO), key_se='FASTQ_SE', drop_empty=True, read_type=params.read_type))
        bowtie2_map.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(bowtie2_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))

        # Initialize tools
        samtools_view = SamtoolsView(Camel.get_instance())
        samtools_sort = SamtoolsSort(Camel.get_instance())
        samtools_sort.update_parameters(threads=threads)
        pipeutils.run_as_pipe([bowtie2_map, samtools_view, samtools_sort], Path(params.dir_))

        # Save output
        SnakemakeUtils.dump_tool_output(samtools_sort, 'BAM', Path(output.BAM))
        bowtie2_map.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_object(bowtie2_map.informs, Path(output.INFORMS))

rule assembly_spades_samtools_depth:
    """
    Runs samtools depth on the BAM file of the reads mapped to the assembly.
    """
    input:
        BAM = rules.assembly_spades_bt2_map.output.BAM
    output:
        TSV = Path(config['working_dir']) / 'assembly_spades' / 'samtools_depth' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_DEPTH_INFORMS
    params:
        dir_ = Path(config['working_dir']) / 'assembly_spades' / 'samtools_depth'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth(camel)
        step = Step(str(rule), samtools_depth, camel, params.dir_)
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        step.run_step()
        samtools_depth.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)
