from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import polish_assembly_short


rule polishing_samtools_index_polypolish:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA = lambda wildcards: Path(config['working_dir']) / str(polish_assembly_short.INPUT_ASSEMBLY_FASTA).format(assembly_type=wildcards.assembly_type)
    output:
        FASTA = Path(config['working_dir']) / polish_assembly_short.OUTPUT_POLISHING_FASTA_INDEX_POLYPOLISH
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'polish' / 'short_reads' / wildcards.assembly_type / 'polypolish'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(samtools, input)
        step = Step(str(rule), samtools, Camel.get_instance(), Path(str(params.dir_)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools, output)

rule polishing_bwa_index:
    """
    Creates a bwa index for the assembly.
    """
    input:
        FASTA_REF = rules.polishing_samtools_index_polypolish.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'polish' / 'short_reads' / '{assembly_type}' / 'polypolish' / 'genome_prefix.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'polish' / 'short_reads' / wildcards.assembly_type / 'polypolish'
    run:
        from camel.app.tools.bwa.bwaindex import BWAIndex
        bwa = BWAIndex(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(bwa, input)
        step = Step(str(rule), bwa, Camel.get_instance(), Path(str(params.dir_)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bwa, output)

rule polishing_read_mapping_1:
    """
    Maps the forward reads against the assembly.
    """
    input:
        FQ_dict = Path(config['working_dir']) / 'fq_dict.io',
        FASTA = rules.polishing_samtools_index_polypolish.output.FASTA,
        INDEX_GENOME_PREFIX_SAMTOOLS = rules.polishing_samtools_index_polypolish.output.FASTA,
        INDEX_GENOME_PREFIX = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'polish' / 'short_reads' / '{assembly_type}' / 'read_mapping' / 'forward' / 'bwa_readmap.io',
        INFORMS = Path(config['working_dir']) / 'polish' / 'short_reads' / '{assembly_type}' / 'read_mapping' / 'forward' / 'commands.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'polish' / 'short_reads' / wildcards.assembly_type / 'read_mapping' / 'forward'
    threads: 8
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.bwa.bwamap import BWAMap
        bwa_map = BWAMap(Camel.get_instance())
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict), 'illumina')
        bwa_map.add_input_files({'FASTQ_SE': [fq_in.pe[0]]})
        bwa_map.update_parameters(threads=threads, all_alns=True)
        SnakemakeUtils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step = Step(str(rule), bwa_map, Camel.get_instance(), Path(str(params.dir_)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bwa_map, output)

rule polishing_read_mapping_2:
    """
    Maps the reverse reads against the assembly.
    """
    input:
        FQ_dict = Path(config['working_dir']) / 'fq_dict.io',
        FASTA = rules.polishing_samtools_index_polypolish.output.FASTA,
        INDEX_GENOME_PREFIX = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'polish' / 'short_reads' / '{assembly_type}' / 'read_mapping' / 'reverse' / 'bwa_readmap.io',
        INFORMS = Path(config['working_dir']) / 'polish' / 'short_reads' / '{assembly_type}' / 'read_mapping' / 'reverse' / 'commands.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'polish' / 'short_reads' / wildcards.assembly_type / 'read_mapping' / 'reverse'
    threads: 8
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.bwa.bwamap import BWAMap
        bwa_map = BWAMap(Camel.get_instance())
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict), 'illumina')
        bwa_map.add_input_files({'FASTQ_SE': [fq_in.pe[1]]})
        bwa_map.update_parameters(threads=threads, all_alns=True)
        SnakemakeUtils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step = Step(str(rule), bwa_map, Camel.get_instance(), Path(str(params.dir_)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bwa_map, output)

rule polishing_polypolish_insert_filter:
    """
    Exclude alignments based on insert size. 
    """
    input:
        SAM_1 = rules.polishing_read_mapping_1.output.SAM,
        SAM_2 = rules.polishing_read_mapping_2.output.SAM
    output:
        SAM = Path(config['working_dir']) / 'polish' / 'short_reads' / '{assembly_type}' / 'read_mapping' / 'alignment_filtered_sam.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'polish' / 'short_reads' / wildcards.assembly_type / 'read_mapping'
    threads: 8
    run:
        from camel.app.tools.polypolish.polypolishinsertfilter import PolypolishInsertFilter
        insert_filter = PolypolishInsertFilter(Camel.get_instance())
        input_sam1 = SnakemakeUtils.load_object(Path(input.SAM_1))
        input_sam2 = SnakemakeUtils.load_object(Path(input.SAM_2))
        insert_filter.add_input_files({'SAM': [input_sam1[0], input_sam2[0]]})
        step = Step(str(rule), insert_filter, Camel.get_instance(), Path(str(params.dir_)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(insert_filter, output)

rule polishing_polypolish:
    """
    First polishing with polypolish.
    """
    input:
        SAM = rules.polishing_polypolish_insert_filter.output.SAM,
        FASTA = rules.polishing_samtools_index_polypolish.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / polish_assembly_short.OUTPUT_POLYPOLISH_FASTA,
        INFORMS = Path(config['working_dir']) / polish_assembly_short.OUTPUT_POLYPOLISH_INFORMS
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'polish' / 'short_reads' / wildcards.assembly_type / 'polypolish',
        polypolish_options = config.get('polishing', {}).get('polypolish', {})
    run:
        from camel.app.tools.polypolish.polypolish import Polypolish
        polypolish = Polypolish(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(polypolish, input)
        polypolish.update_parameters(**params.polypolish_options)
        step = Step(str(rule), polypolish, Camel.get_instance(), Path(str(params.dir_)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(polypolish, output)

rule polishing_samtools_index_polca:
    """
    Creates a samtools index for the polypolish assembly.
    """
    input:
        FASTA = rules.polishing_polypolish.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / polish_assembly_short.OUTPUT_POLISHING_FASTA_INDEX_POLCA
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'polish' / 'short_reads' / wildcards.assembly_type / 'polca'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(samtools, input)
        step = Step(str(rule), samtools, Camel.get_instance(), Path(str(params.dir_)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools, output)

rule polishing_polca:
    """
    Then polishing with Polca.
    """
    input:
        FQ_dict = Path(config['working_dir']) / 'fq_dict.io',
        FASTA = rules.polishing_samtools_index_polca.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / 'polish' / 'short_reads' / '{assembly_type}' / 'polca' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / 'polish' / 'short_reads' / '{assembly_type}' / 'polca' / 'informs.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'polish' / 'short_reads' / wildcards.assembly_type / 'polca',
        polca_options = config.get('polishing', {}).get('polca', {})
    threads: 8
    run:
        from camel.app.tools.polca.polca import Polca
        polca = Polca(Camel.get_instance())
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict),'illumina')
        polca.add_input_files({'FASTQ_PE': fq_in.pe})
        SnakemakeUtils.add_pickle_input(polca, 'FASTA', Path(input.FASTA))
        polca.update_parameters(**params.polca_options, threads=threads)
        step = Step(str(rule), polca, Camel.get_instance(), Path(str(params.dir_)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(polca, output)
