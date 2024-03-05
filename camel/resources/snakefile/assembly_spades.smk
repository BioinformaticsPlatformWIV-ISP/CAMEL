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
        FASTA_Contig = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA,
        INFORMS = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_INFORMS
    params:
        dir_ = Path(config['working_dir']) / 'assembly' / 'spades',
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
