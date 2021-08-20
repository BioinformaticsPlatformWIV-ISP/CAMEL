from pathlib import Path
from camel.app.pipeline.step import Step

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.influenzapipeline.snakefile import assembly
from camel.scripts.influenzapipeline.snakefile import genometyping_blastn

camel = Camel.get_instance()

rule run_assembly:
    """
    Runs the assembly
    """
    input:
        FASTQ = Path(config['working_dir']) / genometyping_blastn.OUTPUT_SEQTK_SUBSAMPLE_FASTQ
    output:
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA
    params:
        running_dir=Path(config['working_dir']) / 'assembly'
    threads: 6
    run:
        from camel.app.components.workflows.assemblywrapper import AssemblyWrapper
        from camel.app.components.workflows.utils.fastqinput import FastqInput

        wrapper = AssemblyWrapper(working_dir=params.running_dir)
        if config['analysis_type'] != 'assembly':
            fq_input = FastqInput.from_fq_dict(input.IO, config['read_type'])
        else:
            fq_files = SnakemakeUtils.load_object(input.FASTQ)
            fq_input = FastqInput(read_type=config['read_type'], pe=fq_files, is_pe=True)
        wrapper.run(config['sample_name'], fq_input, calc_qc_stats=True)
        SnakemakeUtils.dump_object([ToolIOFile(wrapper.output.fasta_contigs)], output.FASTA)

rule create_assembly_blast_db:
    """
    Creates a blast database from the assembly
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA
    output:
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA_BLAST_DB
    params:
        running_dir=Path(config['working_dir']) / 'assembly'
    threads: 6
    run:
        from camel.app.tools.blast.makeblastdb import MakeBlastDb

        makedb = MakeBlastDb(camel)
        SnakemakeUtils.add_pickle_inputs(makedb, input)
        step = Step(str(rule), makedb, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(makedb, output)
