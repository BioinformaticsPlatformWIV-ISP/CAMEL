"""
This Snakefile checks for the presence / absence of the region of difference 1 (RD1) and the csb gene
in Mycobacterium samples.
"""
from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.mycobacteriumpipeline.snakefile import csb_rd

rule rd_csb_report:
    """
    This rule generates the output HTML report for the RD - csb species identification.
    """
    input:
        FASTA = Path(config['working_dir']) / 'gene_detection' / 'csb_rd' / 'db_manager' / 'fasta.io',
        HITS = Path(config['working_dir']) / 'gene_detection' / 'csb_rd' / 'hit_selection' / 'selected-hits.io',
        INFORMS_columns = Path(config['working_dir']) / 'gene_detection' / 'csb_rd' / 'report' / 'informs-columns.io'
    output:
        VAL_HTML = Path(config['working_dir']) / csb_rd.OUTPUT_CSB_RD_REPORT,
        INFORMS = Path(config['working_dir']) / 'csb_rd' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'csb_rd',
        detection_method = config['detection_method']
    run:
        from camel.app.tools.pipelines.mycobacterium.rdcsbreporter import RdCsbReporter
        reporter = RdCsbReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, Camel.get_instance(), Path(params.dir_))
        reporter.update_parameters(hit_type=params.detection_method)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule rd_csb_report_empty:
    """
    This generates an empty RD1 / csb report. 
    """
    output:
        VAL_HTML = Path(config['working_dir']) / csb_rd.OUTPUT_CSB_RD_REPORT_EMPTY
    params:
        dir_ = Path(config['working_dir']) / 'contamination_check' / 'report'
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.mycobacterium.rdcsbreporter import RdCsbReporter
        section = RdCsbReporter.generate_empty_section()
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.VAL_HTML))

rule rd_csb_dump_summary_info:
    """
    Dumps the summary information from the csb / RD1 workflow in tabular format.
    """
    input:
        INFORMS_report = rules.rd_csb_report.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / csb_rd.OUTPUT_CSB_RD_SUMMARY
    params:
        dir_ = Path(config['working_dir']) / 'csb_rd'
    run:
        informs = SnakemakeUtils.load_object(Path(input.INFORMS_report))
        summary_data = [
            ('csb_detected', 'csb' in informs['loci_detected']),
            ('RD1_detected', 'RD1' in informs['loci_detected']),
            ('RD9_detected', 'RD9' in informs['loci_detected']),
            ('rd_csb_species', informs['species'])
        ]
        with open(output.TSV, 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')
