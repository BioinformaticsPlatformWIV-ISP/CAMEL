from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

camel = Camel.get_instance()

rule polishing_copy_fasta:
    """
    Copies the medaka output FASTA file into the polypolish input folder.
    """
    input:
        FASTA = Path(config['working_dir']) / 'medaka' / 'consensus.fasta'
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'input_genome.fasta'
    shell:
        """
        cp {input.FASTA} {output.FASTA}
        """

rule polishing_samtools_index_polypolish:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA_REF = rules.polishing_copy_fasta.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'input_genome.fasta.fai'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polypolish'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex(camel)
        samtools.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule),samtools,camel,params.running_dir,config)
        step.run_step()

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
        step = Step(str(rule), bwa, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bwa, output)

rule polishing_read_mapping:
    """
    Maps the reads against the assembly.
    """
    input:
        FQ_dict = Path(config['working_dir']) / 'trimming' / 'illumina' / 'fq_dict.io',
        INDEX_GENOME_PREFIX_BWA = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX,
        INDEX_GENOME_PREFIX_SAMTOOLS = rules.polishing_samtools_index_polypolish.output.INDEX_GENOME_PREFIX,
        FASTA = rules.polishing_copy_fasta.output.FASTA,
        INDEX_GENOME_PREFIX = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'polishing' / 'read_mapping' / 'bwa_readmap.sam'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'read_mapping'
    threads: 8
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.bwa.bwamap import BWAMap
        bwa_map = BWAMap(camel)
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict), 'illumina')
        bwa_map.add_input_files({'FASTQ_PE': fq_in.pe})
        bwa_map.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step = Step(str(rule), bwa_map, camel, params.running_dir)
        step.run_step()

rule polishing_polypolish:
    """
    First polishing with polypolish.
    """
    input:
        SAM = rules.polishing_read_mapping.output.SAM,
        FASTA = rules.polishing_copy_fasta.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta',
        INFORMS = Path(config['working_dir']) / 'polishing' / 'polypolish'  / 'polypolish.io'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polypolish',
        polypolish_options = config.get('polishing', {}).get('polypolish', {})
    run:
        from camel.app.tools.polypolish.polypolish import Polypolish
        polypolish = Polypolish(camel)
        polypolish.add_input_files({'SAM': [ToolIOFile(Path(input.SAM))], 'FASTA':[ToolIOFile(Path(input.FASTA))]})
        polypolish.update_parameters(**params.polypolish_options)
        step = Step(str(rule), polypolish, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_object(polypolish.informs, Path(output.INFORMS))

rule polishing_copy_fasta_polca:
    input:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta'
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polca' / 'input_genome.fasta'
    shell:
        """
        cp {input.FASTA} {output.FASTA}
        """

rule polishing_samtools_index_polca:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA_REF = rules.polishing_copy_fasta_polca.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'polishing' / 'polca' / 'input_genome.fasta.fai'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polca'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex(camel)
        samtools.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule),samtools,camel,params.running_dir,config)
        step.run_step()

rule polishing_polca:
    """
    Then polishing with Polca.
    """
    input:
        FQ_dict = Path(config['working_dir']) / 'trimming' / 'illumina' / 'fq_dict.io',
        FASTA = rules.polishing_copy_fasta_polca.output.FASTA,
        INDEX = rules.polishing_samtools_index_polca.output.INDEX_GENOME_PREFIX
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polca' / 'input_genome.fasta.PolcaCorrected.fa',
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
        polca.update_parameters(**params.polca_options, threads=threads)
        step = Step(str(rule), polca, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_object(polca.informs, Path(output.INFORMS))

rule polishing_rename_polca_output:
    """
    Renames the fasta file generated by polca.
    """
    input:
        FASTA = rules.polishing_polca.output.FASTA
    output:
         FASTA = Path(config['working_dir']) / 'polishing' / 'polca' / 'polished.fasta'
    shell:
        """
        mv {input.FASTA} {output.FASTA}
        """
