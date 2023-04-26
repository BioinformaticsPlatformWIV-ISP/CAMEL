import pickle
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

camel = Camel.get_instance()

rule copy_fasta:
    input:
        FASTA = Path(config['working_dir']) / 'medaka' / 'consensus.fasta'
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'consensus.fasta'
    shell:
        """
        cp {input.FASTA} {output.FASTA}
        """

rule samtools_index_polypolish:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA_REF = Path(config['working_dir']) / 'medaka' / 'consensus.fasta'
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'medaka' / 'consensus.fasta.fai'
    params:
        running_dir = Path(config['working_dir']) / 'medaka'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex(camel)
        samtools.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule),samtools,camel,params.running_dir,config)
        step.run_step()

rule bwa_index:
    """
    Creates a bwa index for the assembly.
    """
    input:
        FASTA_REF = Path(config['working_dir']) / 'medaka' / 'consensus.fasta'
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'medaka' / 'genome_prefix.io'
    params:
        running_dir = Path(config['working_dir']) / 'medaka'
    run:
        from camel.app.tools.bwa.bwaindex import BWAIndex
        bwa = BWAIndex(camel)
        bwa.add_input_files({'FASTA_REF': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule), bwa, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bwa, output)

rule read_mapping:
    """
    Maps the reads against the assembly.
    """
    input:
        FQ_1P=Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_1P.fastq.gz",
        FQ_2P=Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_2P.fastq.gz",
        INDEX_GENOME_PREFIX_BWA=rules.bwa_index.output.INDEX_GENOME_PREFIX,
        INDEX_GENOME_PREFIX_SAMTOOLS=rules.samtools_index_polypolish.output.INDEX_GENOME_PREFIX,
        FASTA=Path(config['working_dir']) / 'medaka' / 'consensus.fasta',
        INDEX_GENOME_PREFIX = rules.bwa_index.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'polishing' / 'read_mapping' / 'bwa_readmap.sam'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'read_mapping'
    threads: 8
    run:
        from camel.app.tools.bwa.bwamap import BWAMap
        bwa_map = BWAMap(camel)
        bwa_map.add_input_files({'FASTQ_PE': [ToolIOFile(Path(input.FQ_1P)), ToolIOFile(Path(input.FQ_2P))]})
        bwa_map.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step = Step(str(rule), bwa_map, camel, params.running_dir)
        step.run_step()

rule polypolish_polishing:
    """
    First polishing with polypolish.
    """
    input:
        SAM = rules.read_mapping.output.SAM,
        FASTA = Path(config['working_dir']) / 'medaka' / 'consensus.fasta'
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta',
        INFORMS = Path(config['working_dir']) / 'polishing' / 'polypolish'  / 'polypolish.io'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polypolish',
        polypolish_options = config.get('polishing', {}).get('polypolish', {})
    threads: 8
    priority: 1
    run:
        from camel.app.tools.polypolish.polypolish import Polypolish
        polypolish = Polypolish(camel)
        polypolish.add_input_files({'SAM': [ToolIOFile(Path(input.SAM))], 'FASTA':[ToolIOFile(Path(input.FASTA))]})
        polypolish.update_parameters(**params.polypolish_options)
        step = Step(str(rule), polypolish, camel, params.running_dir)
        step.run_step()
        with open(output.INFORMS, 'wb') as handle:
            pickle.dump(polypolish.informs, handle)

rule copy_fasta_polca:
    input:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta'
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polca' / 'consensus.fasta'
    shell:
        """
        cp {input.FASTA} {output.FASTA}
        """

rule samtools_index_polca:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA_REF = rules.copy_fasta_polca.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'polishing' / 'polca' / 'consensus.fasta.fai'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polca'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex(camel)
        samtools.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule),samtools,camel,params.running_dir,config)
        step.run_step()

rule polca_polishing:
    """
    Then polishing with Polca.
    """
    input:
        FQ_1P = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_1P.fastq.gz",
        FQ_2P = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_2P.fastq.gz",
        FASTA = rules.copy_fasta_polca.output.FASTA,
        INDEX = rules.samtools_index_polca.output.INDEX_GENOME_PREFIX
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polca' / 'consensus.fasta.PolcaCorrected.fa',
        INFORMS = Path(config['working_dir']) / 'polishing' / 'polca' / 'polca.io'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polca',
        polca_options = config.get('polishing', {}).get('polca', {})
    threads: 8
    run:
        from camel.app.tools.polca.polca import Polca
        polca = Polca(camel)
        polca.add_input_files({'FASTQ_PE': [ToolIOFile(Path(input.FQ_1P)), ToolIOFile(Path(input.FQ_2P))],
                               'FASTA': [ToolIOFile(Path(input.FASTA))]})
        polca.update_parameters(**params.polca_options)
        polca.update_parameters(threads=threads)
        step = Step(str(rule), polca, camel, params.running_dir, config)
        step.run_step()
        with open(output.INFORMS, 'wb') as handle:
            pickle.dump(polca.informs, handle)

rule samtools_index_short:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA_REF = rules.polca_polishing.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'polishing' / 'polca' / 'consensus.fasta.PolcaCorrected.fa.fai'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polca'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex(camel)
        samtools.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule),samtools,camel,params.running_dir,config)
        step.run_step()
