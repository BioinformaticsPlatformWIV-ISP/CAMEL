from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

from camel.scripts.broadwgs.snakefile import alignment

camel = Camel.get_instance()

rule samtools_convert_to_cram:
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        CRAM = Path(config['working_dir']) / "bamtocram" / "convert" / "cram.io",
    params:
        working_dir = Path(config['working_dir']) / "bamtocram" / "convert",
        output_file = f'{config["sample"]}.cram'
    threads: config['params_smk']['threads']
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        Path(params.working_dir).mkdir(exist_ok=True)

        bam_to_cram = SamtoolsView(camel)
        SnakemakeUtils.add_pickle_input(bam_to_cram,"BAM",input.BAM)
        SnakemakeUtils.add_pickle_input(bam_to_cram,"FASTA_REF",input.FASTA_REF)
        step = Step(rule, bam_to_cram, camel, params.working_dir, config)
        bam_to_cram.update_parameters(
            output_filename = params.output_file,
            output_format = 'CRAM',
            threads = threads,
            **config['rule_params']['bam_to_cram'][rule]
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(bam_to_cram, 'CRAM', output.CRAM)

rule checksum_cram:
    input:
        CRAM = rules.samtools_convert_to_cram.output.CRAM,
    output:
        CRAM_checksum = Path(config['working_dir']) / "bamtocram" / "checksum" / "cram.md5",
    run:
        import subprocess

        cram_file = SnakemakeUtils.load_object(input.CRAM)[0]
        subprocess.run(f"md5sum {cram_file} > {output.CRAM_checksum}", shell = True, executable="/bin/bash")

rule samtools_index_cram:
    input:
        CRAM = rules.samtools_convert_to_cram.output.CRAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        CRAI = Path(config['working_dir']) / "bamtocram" / "index" / "crai.io",
    params:
        working_dir = Path(config['working_dir']) / "bamtocram" / "index",
        output_file = f'{config["sample"]}.crai'
    run:
        from camel.app.tools.samtools.samtoolsindexcram import SamtoolsIndexCram

        Path(params.working_dir).mkdir(exist_ok=True)

        cram_index = SamtoolsIndexCram(camel)
        SnakemakeUtils.add_pickle_inputs(cram_index, input)
        step = Step(rule, cram_index, camel, params.working_dir, config)
        cram_index.update_parameters(output_filename = params.working_dir / params.output_file)
        step.run_step()
        SnakemakeUtils.dump_tool_output(cram_index, 'CRAI', output.CRAI)


rule picard_validate_cram:
    input:
        CRAM = rules.samtools_convert_to_cram.output.CRAM,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        TXT_metrics = Path(config['working_dir']) / "bamtocram" / "metrics" / 'cram_validation_report.io',
    params:
        working_dir = Path(config['working_dir']) / "bamtocram" / "metrics",
    run:
        from camel.app.tools.picard.validatesamfile import ValidateSamFile

        Path(params.working_dir).mkdir(exist_ok=True)

        val_cram = ValidateSamFile(camel)
        SnakemakeUtils.add_pickle_input(val_cram, "BAM", input.CRAM)
        SnakemakeUtils.add_pickle_input(val_cram, "FASTA_REF", input.FASTA_REF)
        step = Step(rule, val_cram, camel, params.working_dir, config)
        val_cram.update_parameters(
            **config['rule_params']['bam_to_cram'][rule]
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(val_cram, 'TXT_report', output.TXT_metrics)
