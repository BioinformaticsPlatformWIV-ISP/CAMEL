from pathlib import Path

from camel.resources.snakefile import gene_detection
from camel.scripts.stecpipeline.snakefile import serotype_detection


rule serotype_detection_run:
    """
    Retrieves the serotype based on the H and O typing. 
    """
    input:
        HITS_O = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_ALL_HITS).format(db='serotype_o'),
        HITS_H = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_ALL_HITS).format(db='serotype_h')
    output:
        VAL_serotype = Path(config['working_dir']) / serotype_detection.OUTPUT_VAL_SEROTYPE
    params:
        running_dir = Path(config['working_dir']) / 'serotype_detection'
    run:
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        from camel.app.tools.pipelines.stec.serotypedetector import SerotypeDetectorEcoli
        from camel.app.pipeline.step import Step
        detector = SerotypeDetectorEcoli(camel)
        SnakemakeUtils.add_pickle_inputs(detector, input)
        step = Step(rule, detector, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(detector, output)

rule serotype_detection_report:
    """
    Creates a simple report section for the detected serotype.
    """
    input:
        VAL_serotype = Path(config['working_dir']) / serotype_detection.OUTPUT_VAL_SEROTYPE if 'serotype' in config['analyses'] else []
    output:
        HTML = Path(config['working_dir']) / serotype_detection.OUTPUT_SEROTYPE_REPORT
    run:
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        from camel.app.components.html.htmlreportsection import HtmlReportSection
        from camel.app.io.tooliovalue import ToolIOValue
        try:
            serotype = SnakemakeUtils.load_object(input.VAL_serotype)[0].value
        except AttributeError:
            serotype = 'NA'
        section = HtmlReportSection(None)
        section.add_paragraph("Detected serotype: <b>{}</b>".format(serotype))
        SnakemakeUtils.dump_object([ToolIOValue(section)], output.HTML)

rule serotype_detection_dump_summary_info:
    """
    Dumps the summary information from the serotype detection.
    """
    input:
        VAL_serotype=Path(config['working_dir']) / serotype_detection.OUTPUT_VAL_SEROTYPE if 'serotype' in config['analyses'] else []
    output:
        Path(config['working_dir']) / serotype_detection.OUTPUT_SEROTYPE_SUMMARY
    run:
        serotype = SnakemakeUtils.load_object(input.VAL_serotype)[0].value
        with open(output[0], 'w') as handle:
            handle.write(f'serotype\t{serotype}')
            handle.write('\n')
