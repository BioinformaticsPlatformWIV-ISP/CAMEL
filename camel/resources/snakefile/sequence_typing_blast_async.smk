from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.assembly_spades import OUTPUT_ASSEMBLY_FASTA

camel = Camel.get_instance()

rule typing_blast_async:
    """
    Performs asynchronous sequence typing using BLAST+.
    """
    input:
        FASTA = Path(config['working_dir']) / OUTPUT_ASSEMBLY_FASTA,
        DIR = lambda wildcards: config['sequence_typing'][wildcards.scheme]['path']
    output:
        IO = Path(config['working_dir']) / 'typing' / '{scheme}' / '{locus_type}' / '{detection_method}' / 'all-hits.io'
    params:
        dir_working = lambda wildcards: Path(config['working_dir']) / 'typing' / wildcards.scheme / wildcards.locus_type
    threads: 16
    run:
        from camel.app.tools.pipelines.sequence_typing.typeasync import TypeAsync

        # Create working directory
        dir_working = Path(str(params.dir_working))
        dir_working.mkdir(parents=True, exist_ok=True)

        # Run the tool
        typer = TypeAsync(Camel.get_instance())
        typer.add_input_files({
            'FASTA': SnakemakeUtils.load_object(Path(input.FASTA)),
            'DIR': [ToolIODirectory(Path(str(input.DIR)))]
        })
        typer.update_parameters(threads=threads)
        typer.run(dir_working)

        # Save the tool output
        SnakemakeUtils.dump_object(typer.tool_outputs['VAL_hits'], Path(output.IO))
