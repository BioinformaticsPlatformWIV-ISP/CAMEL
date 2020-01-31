"""
This Snakefile is used to perform a contamination check using KRAKEN. A report is generated with statistics and an
interactive KRONA visualization.
"""
from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import contamination_check_kraken


camel = Camel.get_instance()


rule contamination_check_get_db:
    """
    Converts the Kraken database to a pickle IO object.
    """
    output:
        DB = Path(config['working_dir']) / 'contamination_check' / 'db.io'
    params:
        db = config['contamination_check']['db']
    run:
        from camel.app.io.tooliodirectory import ToolIODirectory
        SnakemakeUtils.dump_object([ToolIODirectory(params.db)], output.DB)

rule contamination_check_kraken_run:
    """
    Assign taxonomic labels to reads using KRAKEN.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
        DB = rules.contamination_check_get_db.output.DB
    output:
        TSV = Path(config['working_dir']) / 'contamination_check' / 'kraken' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'contamination_check' / 'kraken'
    threads: 8
    priority: 1
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.kraken.kraken import Kraken
        kraken = Kraken(camel)
        SnakemakeUtils.add_pickle_inputs(kraken, input, ['DB'])
        kraken.add_input_files(SnakePipelineUtils.extracts_fq_input(input.IO, key_se='FASTQ', drop_se=True))
        step = Step(rule, kraken, camel, params.running_dir, config)
        kraken.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(kraken, output)

rule contamination_check_kraken_report:
    """
    Creteas a report based on the Kraken output.
    """
    input:
        TSV = rules.contamination_check_kraken_run.output.TSV,
        DB = rules.contamination_check_get_db.output.DB
    output:
        TSV = Path(config['working_dir']) / 'contamination_check' / 'kraken' / 'tsv-report.io'
    params:
        running_dir = Path(config['working_dir']) / 'contamination_check' / 'kraken'
    run:
        from camel.app.tools.kraken.krakenreport import KrakenReport
        kraken_report = KrakenReport(camel)
        SnakemakeUtils.add_pickle_inputs(kraken_report, input)
        step = Step(rule, kraken_report, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(kraken_report, output)

rule contamination_check_kraken_report_parser:
    """
    Parses the Kraken report and looks for contamination at the species level. 
    """
    input:
        TSV = rules.contamination_check_kraken_report.output.TSV
    output:
        INFORMS = Path(config['working_dir']) / 'contamination_check' / 'kraken' / 'informs-contamination.io'
    params:
        running_dir = Path(config['working_dir']) / 'contamination_check' / 'kraken',
        expected_species = config['contamination_check']['expected_species']
    run:
        from camel.app.tools.kraken.krakenreportparser import KrakenReportParser
        report_parser = KrakenReportParser(camel)
        SnakemakeUtils.add_pickle_inputs(report_parser, input)
        step = Step(rule, report_parser, camel, params.running_dir, config)
        report_parser.update_parameters(expected_species=params.expected_species)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(report_parser, output)

rule contamination_check_krona:
    """
    Creates an interactive pie chart displaying the Kraken output.
    """
    input:
        TSV = rules.contamination_check_kraken_run.output.TSV
    output:
        HTML = Path(config['working_dir']) / 'contamination_check' / 'krona' / 'html.io'
    params:
        running_dir = Path(config['working_dir']) / 'contamination_check' / 'krona'
    run:
        from camel.app.tools.krona.krona import Krona
        krona = Krona(camel)
        SnakemakeUtils.add_pickle_inputs(krona, input)
        step = Step(rule, krona, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(krona, output)

rule contamination_check_report:
    """
    Creates a report containing the results of the contamination check.
    """
    input:
        HTML_Krona = rules.contamination_check_krona.output.HTML,
        INFORMS_species = rules.contamination_check_kraken_report_parser.output.INFORMS,
        INFORMS_kraken = rules.contamination_check_kraken_run.output.INFORMS,
        TSV = rules.contamination_check_kraken_report.output.TSV,
        DB = rules.contamination_check_get_db.output.DB
    output:
        VAL_HTML = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'contamination_check' / 'report'
    run:
        from camel.app.tools.pipelines.quality_checks.htmlreportercontamination import HtmlReporterContamination
        reporter = HtmlReporterContamination(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule contamination_check_report_empty:
    """
    Generates an empty contamination check report.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY
    params:
        running_dir = Path(config['working_dir']) / 'contamination_check' / 'report'
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.quality_checks.htmlreportercontamination import HtmlReporterContamination
        section = HtmlReporterContamination.generate_empty_section()
        SnakemakeUtils.dump_object([ToolIOValue(section)], output.VAL_HTML)

rule contamination_check_dump_summary_info:
    """
    Dumps the summary information for the contamination check in tabular format.
    """
    input:
        INFORMS_species = rules.contamination_check_kraken_report_parser.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_SUMMARY
    run:
        informs = SnakemakeUtils.load_object(input.INFORMS_species)
        summary_data = [
            ('expected_species', informs['expected'][0]),
            ('expected_species_occurrence', informs['expected'][1]),
            ('contaminants_warn', str(informs['contaminants_warn'])),
            ('contaminants_fail', str(informs['contaminants_fail']))
        ]
        with open(output[0], 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')
