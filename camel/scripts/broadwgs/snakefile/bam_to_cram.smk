from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue

camel = Camel()

rule samtools_convert_to_cram:
    input:
        BAM = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam" / (config["sample"] + ".bam.io"),
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        CRAM = Path(config['working_dir']) / "bamtocram" / "convert" / (config["sample"] + ".cram.io"),
    params:
        working_dir = Path(config['working_dir']) / "bamtocram" / "convert"
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        Path(params.working_dir).mkdir(exist_ok=True)

        bamtocram = SamtoolsView(camel)
        SnakemakeUtils.add_pickle_input(bamtocram,"BAM",input.BAM)
        step = Step(rule, bamtocram, camel, params.working_dir, config)
        bamtocram.update_parameters(
            output_filename = "samtools_view.cram",
            cram_out = ""
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(bamtocram, 'BAM', output.CRAM)

rule checksum_cram:
    input:
        CRAM = Path(config['working_dir']) / "bamtocram" / "convert" / (config["sample"] + ".cram.io"),
    output:
        CRAM_checksum = Path(config['working_dir']) / "bamtocram" / (config["sample"] + ".cram.md5"),
    run:
        import subprocess

        cram_file = SnakemakeUtils.load_object(input.CRAM)[0]
        subprocess.run(f"md5sum {cram_file} > {output.CRAM_checksum}", shell = True, executable="/bin/bash")

rule samtools_index_cram:
    input:
        CRAM = Path(config['working_dir']) / "bamtocram" / "convert" / (config["sample"] + ".cram.io"),
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        CRAI = Path(config['working_dir']) / "bamtocram" / "index" / (config["sample"] + ".sample.crai.io"),
    params:
        working_dir = Path(config['working_dir']) / "bamtocram" / "index"
    run:
        from camel.app.tools.samtools.samtoolsindexcram import SamtoolsIndexCram

        Path(params.working_dir).mkdir(exist_ok=True)

        cram_index = SamtoolsIndexCram(camel)
        SnakemakeUtils.add_pickle_inputs(cram_index, input)
        step = Step(rule, cram_index, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(cram_index, 'CRAI', output.CRAI)


rule validate_cram:
    input:
        CRAM = Path(config['working_dir']) / "bamtocram" / "convert" / (config["sample"] + ".cram.io"),
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        TXT_metrics = Path(config['working_dir']) / "bamtocram" / "metrics" / (config["sample"] + "_cram_validation_report.io"),
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / "bamtocram" / "metrics",
        max_output = 1000000000,
    run:
        from camel.app.tools.picard.validatesamfile import ValidateSamFile

        Path(params.working_dir).mkdir(exist_ok=True)

        val_cram = ValidateSamFile(camel)
        SnakemakeUtils.add_pickle_input(val_cram, "BAM", input.CRAM)
        SnakemakeUtils.add_pickle_input(val_cram, "FASTA_REF", input.FASTA_REF)
        step = Step(rule, val_cram, camel, params.working_dir, config)
        val_cram.update_parameters(
            max_output = params.max_output,
            ignore = "MISSING_TAG_NM",
            mode = "VERBOSE",
            is_bisulfite_sequenced = "false"
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(val_cram, 'TXT_report', output.TXT_metrics)
