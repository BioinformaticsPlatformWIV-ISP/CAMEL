"""
This Snakefile performs species/subspecies identification of Mycobacterium strains based on the hsp65 gene sequences.
"""
from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.mycobacteriumpipeline.snakefile import hsp65

rule hps65_report:
    """
    Creates a HTML report for the HSP65 (sub)species identification.
    """
    input:
        INFORMS_hits = Path(config['working_dir']) / 'gene_detection' / 'hsp65' / 'hit_selection' / 'selected-hits.io',
        INFORMS_columns = Path(config['working_dir']) / 'gene_detection' / 'hsp65' / 'report' / 'informs-columns.io',
        FASTA_DB = Path(config['working_dir']) / 'gene_detection' / 'hsp65' / 'db_manager' / 'fasta.io'
    output:
        VAL_HTML = Path(config['working_dir']) / hsp65.OUTPUT_HSP65_REPORT,
        INFORMS = Path(config['working_dir']) / 'hsp65' /  'report' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'hsp65',
        hit_type = config['detection_method']
    run:
        from camel.app.tools.pipelines.mycobacterium.hsp65reporter import Hsp65Reporter
        reporter = Hsp65Reporter(Camel.get_instance())
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        reporter.update_parameters(hit_type=params.hit_type)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule hps65_report_empty:
    """
    Creates an empty HTML report for the HSP65 (sub)species identification.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / hsp65.OUTPUT_HSP65_REPORT_EMPTY
    params:
        dir_ = Path(config['working_dir']) / 'contamination_check' / 'report'
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.mycobacterium.hsp65reporter import Hsp65Reporter
        section = Hsp65Reporter.generate_empty_section()
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.VAL_HTML))

rule hsp65_dump_summary_info:
    """
    Dumps the summary information for the hsp65 workflow in tabular format.
    """
    input:
        INFORMS = rules.hps65_report.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / hsp65.OUTPUT_HSP65_SUMMARY
    run:
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        summary_data = [('hsp65_species', informs['hits']),]
        with open(output.TSV, 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')
