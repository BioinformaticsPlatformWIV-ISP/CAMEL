from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.influenzapipeline.snakefile import assembly

camel = Camel.get_instance()

rule run_assembly:
    """
    Runs the assembly
    """
    input:
        IO=Path(config['working_dir']) / 'fq_dict.io'
    output:
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA
    params:
        running_dir=Path(config['working_dir']) / 'assembly'
    threads: 6
    run:
        from camel.app.components.workflows.assemblywrapper import AssemblyWrapper
        from camel.app.components.workflows.utils.fastqinput import FastqInput

        wrapper = AssemblyWrapper(working_dir=params.running_dir)
        fq_input = FastqInput.from_fq_dict(input.IO, config['read_type'])
        wrapper.run(config['sample_name'], fq_input, calc_qc_stats=True)
        SnakemakeUtils.dump_object([ToolIOFile(wrapper.output.fasta_contigs)], output.FASTA)
