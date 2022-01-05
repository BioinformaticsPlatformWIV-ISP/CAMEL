from pathlib import Path

from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.error.snakemakeexecutionerror import SnakemakeExecutionError
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.broadwgs.snakefile import alignment
from camel.scripts.broadwgs import references

camel = Camel.get_instance()

rule samtools_convert_to_cram:
    """
    Convert the aligned BAM to CRAM for storage.
    """
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
        FASTA_REF = Path(config['working_dir']) / references.FASTA_GENOME_FILE,
    output:
        CRAM = Path(config['working_dir']) / "bamtocram" / "convert" / "cram.io",
    params:
        working_dir = Path(config['working_dir']) / "bamtocram" / "convert",
        output_file = f'{config["sample"]}.cram'
    threads: config["params_smk"]["threads_cram"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_cram"]
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        Path(params.working_dir).mkdir(exist_ok=True)

        bam_to_cram = SamtoolsView(camel)
        SnakemakeUtils.add_pickle_input(bam_to_cram, "BAM", Path(input.BAM))
        SnakemakeUtils.add_pickle_input(bam_to_cram, "FASTA_REF", Path(input.FASTA_REF))
        step = Step(rule, bam_to_cram, camel, params.working_dir)
        bam_to_cram.update_parameters(
            output_filename = params.output_file,
            output_format = 'CRAM',
            threads = threads,
            **config['rule_params']['bam_to_cram'][rule]
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(bam_to_cram, 'CRAM', Path(output.CRAM))

rule checksum_cram:
    """
    Calculate the checksum of the CRAM file.
    """
    input:
        CRAM = rules.samtools_convert_to_cram.output.CRAM,
    output:
        CRAM_checksum = Path(config['working_dir']) / "bamtocram" / "convert" / "cram.md5",
    threads: config["params_smk"]["threads_cram"]
    params:
        working_dir = Path(config['working_dir']) / "bamtocram" / "convert"
    resources:
        mem_mb=config["params_smk"]["memory_mb_cram"]
    run:
        cram_file = SnakemakeUtils.load_object(Path(input.CRAM))[0]

        cmd_checksum = Command(f"md5sum {cram_file} > {output.CRAM_checksum}")
        cmd_checksum.run(params.working_dir)
        if cmd_checksum.returncode != 0:
            raise SnakemakeExecutionError(cmd_checksum.stdout, cmd_checksum.stderr)

rule samtools_index_cram:
    """
    Index the CRAM file.
    """
    input:
        CRAM = rules.samtools_convert_to_cram.output.CRAM,
        FASTA_REF = Path(config['working_dir']) / references.FASTA_GENOME_FILE
    output:
        CRAI = Path(config['working_dir']) / "bamtocram" / "convert" / "cram.crai.io",
    params:
        working_dir = Path(config['working_dir']) / "bamtocram" / "convert",
        output_file = f'{config["sample"]}.cram'
    threads: config["params_smk"]["threads_cram"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_cram"]
    run:
        from camel.app.tools.samtools.samtoolsindexcram import SamtoolsIndexCram

        Path(params.working_dir).mkdir(exist_ok=True)

        cram_index = SamtoolsIndexCram(camel)
        SnakemakeUtils.add_pickle_inputs(cram_index, input)
        step = Step(rule, cram_index, camel, params.working_dir)
        cram_index.update_parameters(output_filename = params.working_dir / params.output_file)
        step.run_step()
        SnakemakeUtils.dump_tool_output(cram_index, 'CRAI', Path(output.CRAI))
