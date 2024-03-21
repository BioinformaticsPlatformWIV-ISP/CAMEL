from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

camel = Camel.get_instance()


rule assembly_flye_run:
    """
    De-novo assembly using Flye.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        FASTA = Path(config['working_dir']) / 'assembly' / 'flye' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / 'assembly' / 'flye' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'assembly' / 'flye',
        flye_options = config.get('assembly', {}).get('flye', {})
    threads: 16
    priority: 1
    run:
        from camel.app.tools.flye.flye import Flye
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        flye = Flye(camel)

        # Reformat FASTQ dictionary
        fq_dict = SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_se='FASTQ', read_type='SE')
        flye.add_input_files(fq_dict)
        step = Step(str(rule), flye, camel, params.dir_)
        flye.update_parameters(**params.flye_options)
        flye.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(flye, output)
