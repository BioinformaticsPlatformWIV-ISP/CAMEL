from pathlib import Path
import glob

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue

camel = Camel()

scatter_files = [("temp_" + "{:04d}" + "_of_" + str(config["parameters"]["haplotype_scatter_count"])).format(s) for s in range(1,config["parameters"]["haplotype_scatter_count"]+1)]

rule scatter_interval_list:
    input:
        TXT_intervals = Path(config['working_dir']) / "ref_input" / "calling_intervals.io"
    output:
        TXT_intervalLists = Path(config['working_dir']) / "variant_calling" / "varcalling_intervals" / "interval_lists.io"
    params:
        working_dir = Path(config['working_dir']) / "variant_calling" / "varcalling_intervals",
        break_bands_at_multiples_of = config["parameters"]["break_bands_at_multiples_of"],
        scatter_count = config["parameters"]["haplotype_scatter_count"]
    run:
        from camel.app.tools.picard.intervallisttools import IntervalListTools

        Path(params.working_dir).mkdir(exist_ok=True)

        generate_interval_list = IntervalListTools(camel)
        SnakemakeUtils.add_pickle_inputs(generate_interval_list, input)
        step = Step(rule, generate_interval_list, camel, params.working_dir, config)
        generate_interval_list.update_parameters(
            scatter_count = params.scatter_count,
            subdivision_mode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW",
            unique = "true",
            sort = "true",
            break_bands_at_multiples_of = params.break_bands_at_multiples_of,
            output_dir = params.working_dir
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(generate_interval_list, 'TXT_intervalLists', output.TXT_intervalLists)

rule create_io_intervalLists:
    input:
        TXT_intervalLists = Path(config['working_dir']) / "variant_calling" / "varcalling_intervals" / "interval_lists.io"
    output:
        TXT_intervalList = Path(config['working_dir']) / "variant_calling" / "varcalling_intervals" / "{scatter}" / "scattered.interval_list.io",
    run:
        for interval_list in SnakemakeUtils.load_object(input.TXT_intervalLists):
            intervalList_io = interval_list.path + ".io"
            SnakemakeUtils.dump_object([interval_list], intervalList_io)

rule gatk4_var_caller:
    input:
        BAM = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam" / (config["sample"] + ".bam.io"),
        TXT_intervals = Path(config['working_dir']) / "variant_calling" / "varcalling_intervals" / "{scatter}" / "scattered.interval_list.io",
        FASTA_REF = Path(config['working_dir']) / "ref_input" / "fasta_reference_human_value_file.io",
    output:
        VCF = Path(config['working_dir']) / "variant_calling" / "vcf" / config["sample"] / "{scatter}.vcf.gz.io",
        bamout = Path(config['working_dir']) / "variant_calling" / "vcf" / config["sample"] / "{scatter}.bamout.bam"
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / "variant_calling" / "vcf" / config["sample"] / wildcards.scatter,
        hc_divisor = config["parameters"]["hc_divisor"]
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
            contamination = "0",
            #new_qual = "", ##TODO deprecated
            gqb = "10,20,30,40,50,60,70,80,90",
            annotation_group = "StandardAnnotation,StandardHCAnnotation,AS_StandardAnnotation",
            emit_ref_confidence = "GVCF"
        )
        step.run_step()
        SnakemakeUtils.dump_tool_output(hc, 'VCF', output.VCF)

rule picard_merge_vcfs:
    input:
        VCF = expand(Path(config['working_dir']) / "variant_calling" / "vcf" / config["sample"] / "{scatter}.vcf.gz.io", scatter = scatter_files)
    output:
        VCF = Path(config['working_dir']) / "variant_calling" / "merge_vcf" / (config["sample"] + ".g.vcf.gz.io"),
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / "variant_calling" / "merge_vcf"
    run:
        from camel.app.tools.picard.mergevcfs import MergeVCFs

        merge_vcf = MergeVCFs(camel)
        merge_vcf.add_input_files({"VCF": [SnakemakeUtils.load_object(path)[0] for path in input.VCF]})
        step = Step(rule, merge_vcf, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(merge_vcf, "VCF", output.VCF)



