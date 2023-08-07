import shutil
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import short_read_polishing

camel = Camel.get_instance()


rule polishing_copy_fasta:
    """
    Copies the input FASTA file into the polypolish input folder.
    """
    input:
        FASTA = Path(config['working_dir']) / short_read_polishing.INPUT_ASSEMBLY_FASTA
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'input_genome.fasta'
    run:
        fasta_file = SnakemakeUtils.load_object(Path(str(input.FASTA)))[0].path
        shutil.copyfile(fasta_file, output.FASTA)

rule polishing_samtools_index_polypolish:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA_REF = rules.polishing_copy_fasta.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'samtools_index.io'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polypolish'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex(camel)
        samtools.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule), samtools, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools, output)

rule polishing_bwa_index:
    """
    Creates a bwa index for the assembly.
    """
    input:
        FASTA_REF = rules.polishing_copy_fasta.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'genome_prefix.io'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polypolish'
    run:
        from camel.app.tools.bwa.bwaindex import BWAIndex
        bwa = BWAIndex(camel)
        bwa.add_input_files({'FASTA_REF': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule), bwa, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bwa, output)

rule polishing_read_mapping_1:
    """
    Maps the forward reads against the assembly.
    """
    input:
        FQ_dict = Path(config['working_dir']) / short_read_polishing.INPUT_READS_FASTQ,
        INDEX_GENOME_PREFIX_BWA = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX,
        INDEX_GENOME_PREFIX_SAMTOOLS = rules.polishing_samtools_index_polypolish.output.FASTA,
        FASTA = rules.polishing_copy_fasta.output.FASTA,
        INDEX_GENOME_PREFIX = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'polishing' / 'read_mapping' / 'forward' / 'bwa_readmap.io'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'read_mapping' / 'forward'
    threads: 8
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.bwa.bwamap import BWAMap
        bwa_map = BWAMap(camel)
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict), 'illumina')
        bwa_map.add_input_files({'FASTQ_SE': [fq_in.pe[0]]})
        bwa_map.update_parameters(threads=threads, all_alns=True)
        SnakemakeUtils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step = Step(str(rule), bwa_map, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bwa_map, output)

rule polishing_read_mapping_2:
    """
    Maps the forward reads against the assembly.
    """
    input:
        FQ_dict = Path(config['working_dir']) / short_read_polishing.INPUT_READS_FASTQ,
        INDEX_GENOME_PREFIX_BWA = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX,
        INDEX_GENOME_PREFIX_SAMTOOLS = rules.polishing_samtools_index_polypolish.output.FASTA,
        FASTA = rules.polishing_copy_fasta.output.FASTA,
        INDEX_GENOME_PREFIX = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'polishing' / 'read_mapping' / 'reverse' / 'bwa_readmap.io'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'read_mapping' / 'reverse'
    threads: 8
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.bwa.bwamap import BWAMap
        bwa_map = BWAMap(camel)
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict), 'illumina')
        bwa_map.add_input_files({'FASTQ_SE': [fq_in.pe[1]]})
        bwa_map.update_parameters(threads=threads, all_alns=True)
        SnakemakeUtils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step = Step(str(rule), bwa_map, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bwa_map, output)

rule polishing_polypolish_insert_filter:
    input:
        SAM_1 = rules.polishing_read_mapping_1.output.SAM,
        SAM_2 = rules.polishing_read_mapping_2.output.SAM
    output:
        SAM = Path(config['working_dir']) / 'polishing' / 'read_mapping' / 'alignment_filtered.io'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'read_mapping'
    threads: 8
    run:
        from camel.app.tools.polypolish.polypolishinsertfilter import PolypolishInsertFilter
        insert_filter = PolypolishInsertFilter(camel)
        test_1 = SnakemakeUtils.load_object(Path(input.SAM_1))
        insert_filter.add_input_files({'SAM': [SnakemakeUtils.load_object(Path(input.SAM_1))[0],
                                               SnakemakeUtils.load_object(Path(input.SAM_2))[0]]})
        step = Step(str(rule), insert_filter, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(insert_filter, output)

rule polishing_polypolish:
    """
    First polishing with polypolish.
    """
    input:
        SAM_filtered = rules.polishing_polypolish_insert_filter.output.SAM,
        FASTA = rules.polishing_copy_fasta.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.io',
        INFORMS = Path(config['working_dir']) / 'polishing' / 'polypolish'  / 'polypolish.io'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polypolish',
        polypolish_options = config.get('polishing', {}).get('polypolish', {})
    run:
        from camel.app.tools.polypolish.polypolish import Polypolish
        polypolish = Polypolish(camel)
        polypolish.add_input_files({
            'SAM': SnakemakeUtils.load_object(Path(input.SAM_filtered)),
            'FASTA':[ToolIOFile(Path(input.FASTA))]})
        step = Step(str(rule), polypolish, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(polypolish, output)

rule polishing_copy_fasta_polca:
    input:
        FASTA = rules.polishing_polypolish.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polca' / 'input_genome.fasta'
    run:
        fasta_file = SnakemakeUtils.load_object(Path(str(input.FASTA)))[0].path
        shutil.copyfile(fasta_file, output.FASTA)

rule polishing_samtools_index_polca:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA_REF = rules.polishing_copy_fasta_polca.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polca' / 'samtools_index.io'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polca'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex(camel)
        samtools.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule), samtools, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools, output)

rule polishing_polca:
    """
    Then polishing with Polca.
    """
    input:
        FQ_dict = Path(config['working_dir']) / short_read_polishing.INPUT_READS_FASTQ,
        FASTA = rules.polishing_copy_fasta_polca.output.FASTA,
        INDEX = rules.polishing_samtools_index_polca.output.FASTA
    output:
        # FASTA = Path(config['working_dir']) / 'polishing' / 'polca' / 'input_genome.fasta.PolcaCorrected.fa',
        FASTA = Path(config['working_dir']) / 'polishing' / 'polca' / 'polished.io',
        INFORMS = Path(config['working_dir']) / 'polishing' / 'polca' / 'polca.io'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polca',
        polca_options = config.get('polishing', {}).get('polca', {})
    threads: 8
    run:
        from camel.app.tools.polca.polca import Polca
        polca = Polca(camel)
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict),'illumina')
        polca.add_input_files({'FASTQ_PE': fq_in.pe, 'FASTA': [ToolIOFile(Path(input.FASTA))]})
        polca.update_parameters(threads=threads)
        step = Step(str(rule), polca, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(polca, output)

rule polishing_rename_polca_output:
    """
    Renames the fasta file generated by polca.
    """
    input:
        FASTA = rules.polishing_polca.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polca' / 'polished.fasta'
    run:
        shutil.copyfile(SnakemakeUtils.load_object(Path(str(input.FASTA)))[0].path, output.FASTA)

rule polishing_dump_polca_output:
    """
    Dumps the output fasta file.
    """
    input:
        FASTA = rules.polishing_rename_polca_output.output.FASTA
    output:
         FASTA = Path(config['working_dir']) / short_read_polishing.OUTPUT_POLISHING_FASTA
    run:
        SnakemakeUtils.dump_object(ToolIOFile(Path(input.FASTA)), Path(output.FASTA))
