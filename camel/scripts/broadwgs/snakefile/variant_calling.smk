from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

from camel.scripts.broadwgs.snakefile import alignment

camel = Camel.get_instance()

rule picard_create_interval_lists:
    input:
        TXT_intervals = Path(config['working_dir']) / "ref_input" / "calling_intervals.io"
    output:
        TXT_intervalList = Path(config['working_dir']) / "variant_calling" / "varcalling_intervals" / f'temp_{{scatter}}_of_{config["rule_params"]["variant_calling"]["picard_create_interval_lists"]["scatter_count"]}' / "scattered.interval_list.io"
    params:
        working_dir = Path(config['working_dir']) / "variant_calling" / "varcalling_intervals",
    run:
        from camel.app.tools.picard.intervallisttools import IntervalListTools

        Path(params.working_dir).mkdir(exist_ok=True)

        generate_interval_list = IntervalListTools(camel)
        SnakemakeUtils.add_pickle_inputs(generate_interval_list, input)
        step = Step(rule, generate_interval_list, camel, params.working_dir, config)
        generate_interval_list.update_parameters(
            **config['rule_params']['variant_calling'][rule],
            output = params.working_dir
        )
        step.run_step()

        for interval_list in step.outputs['TXT_intervalLists']:
            intervalList_io = interval_list.path + ".io"
            SnakemakeUtils.dump_object([interval_list], intervalList_io)

rule gatk4_haplotype_caller:
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
        TXT_intervals = rules.picard_create_interval_lists.output.TXT_intervalList,
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        VCF = Path(config['working_dir']) / "variant_calling" / "haplotype_caller" / f'temp_{{scatter}}_of_{config["rule_params"]["variant_calling"]["picard_create_interval_lists"]["scatter_count"]}.vcf.gz.io',
        bamout = Path(config['working_dir']) / "variant_calling" / "haplotype_caller" / f'temp_{{scatter}}_of_{config["rule_params"]["variant_calling"]["picard_create_interval_lists"]["scatter_count"]}.bamout.bam'
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / "variant_calling" / "haplotype_caller" / f'temp_{wildcards.scatter}_of_{config["rule_params"]["variant_calling"]["picard_create_interval_lists"]["scatter_count"]}',
    run:
        import subprocess
        from camel.app.tools.gatk4.gatk4haplotypecaller import GATK4HaplotypeCaller

        Path(params.working_dir).mkdir(exist_ok=True)

        # BAM output not required, but will crash if file not present
        subprocess.run(f"touch {output.bamout}", shell = True, executable="/bin/bash")

        hc = GATK4HaplotypeCaller(camel)
        SnakemakeUtils.add_pickle_inputs(hc, input)
        step = Step(rule, hc, camel, params.working_dir, config)
        hc.update_parameters(
            **config['rule_params']['variant_calling'][rule],
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(hc, 'VCF', output.VCF)

rule picard_merge_vcfs:
    input:
        VCF = expand(rules.gatk4_haplotype_caller.output.VCF, scatter = ["{:04d}".format(s) for s in range(1, (config["rule_params"]["variant_calling"]["picard_create_interval_lists"]["scatter_count"] + 1))])
    output:
        VCF = Path(config['working_dir']) / "variant_calling" / "merge_vcf" / "output.g.vcf.gz.io",
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / "variant_calling" / "merge_vcf"
    run:
        from camel.app.tools.picard.mergevcfs import MergeVCFs

        merge_vcf = MergeVCFs(camel)
        merge_vcf.add_input_files({"VCF": [SnakemakeUtils.load_object(path)[0] for path in input.VCF]})
        step = Step(rule, merge_vcf, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(merge_vcf, "VCF", output.VCF)
