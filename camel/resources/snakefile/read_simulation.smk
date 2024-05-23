from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import assembly

rule read_simulation_fastq_from_fasta:
    """
    Simulates illumina reads from the input FASTA file.
    This rule is executed when the input type is FASTA and no VCF input file has been provided.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.get_fasta_raw(config)
    output:
        FASTQ_PE = Path(config['working_dir']) / 'read_simulation' / 'art' / 'fastq.io',
        INFORMS = Path(config['working_dir']) / 'read_simulation' / 'art' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'read_simulation' / 'art'
    run:
        from camel.app.tools.art.art import ART
        art = ART(Camel.get_instance())
        step = Step(str(rule), art, Camel.get_instance(), params.running_dir)
        SnakemakeUtils.add_pickle_inputs(art, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(art, output)
