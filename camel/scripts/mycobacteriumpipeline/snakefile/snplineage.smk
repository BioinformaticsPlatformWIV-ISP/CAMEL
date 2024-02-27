from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import variant_calling, variant_filtering
from camel.scripts.mycobacteriumpipeline.snakefile import snplineage


rule snp_lineage_detection:
    """
    Detects the SNP lineage based on the variants detected in the sample.
    """
    input:
        VCF = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_UNFILTERED_VCF,
        VCF_filt = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_VCF
    output:
        INFORMS = Path(config['working_dir']) / snplineage.OUTPUT_SNP_LINEAGE_INFORMS
    params:
        dir_ = Path(config['working_dir'], 'snp_lineage'),
        lineage_snp_positions = config['snp_lineage']['bed']
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.tools.pipelines.mycobacterium.snplineagedetector import SNPLineageDetector
        detector = SNPLineageDetector(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(detector, input)
        detector.add_input_files({'BED': [ToolIOFile(Path(params.lineage_snp_positions))]})
        step = Step(str(rule), detector, Camel.get_instance(), Path(params.dir_))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(detector, output)

rule snp_lineage_report:
    """
    Creates the HTML output report for the SNP lineage report.
    """
    input:
        INFORMS_detection = rules.snp_lineage_detection.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / snplineage.OUTPUT_SNP_LINEAGE_REPORT
    params:
        dir_ = Path(config['working_dir'],'snp_lineage'),
        lineage_snp_positions = config['snp_lineage']['bed']
    run:
        from camel.app.tools.pipelines.mycobacterium.snplineagereporter import SNPLineageReporter
        lineage_snp = SNPLineageReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(lineage_snp, input)
        step = Step(str(rule), lineage_snp, Camel.get_instance(), Path(params.dir_))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(lineage_snp, output)

rule snp_lineage_report_empty:
    """
    Creates an empty report when the assay is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / snplineage.OUTPUT_SNP_LINEAGE_REPORT_EMPTY
    run:
        from camel.app.tools.pipelines.mycobacterium.snplineagereporter import SNPLineageReporter
        from camel.app.snakemake.snakepipelineutils import  SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section(SNPLineageReporter.TITLE, Path(output.VAL_HTML), 2)

rule snp_lineage_dump_summary_info:
    """
    Dumps the summary information from the SNP lineage assay.
    """
    input:
        INFORMS = rules.snp_lineage_detection.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / snplineage.OUTPUT_SNP_LINEAGE_SUMMARY
    run:
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        all_lineages = ', '.join([
            d['lineage'].id_ for _, d in informs['detected_lineage_by_level'].items() if d is not None])
        with open(output.TSV, 'w') as handle:
            handle.write(f"snp_lineages\t{all_lineages}")
            handle.write('\n')
