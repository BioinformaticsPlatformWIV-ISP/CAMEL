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


rule contamination_check_kraken2_run:
    """
    Assigns taxonomic labels to reads using KRAKEN2.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
        DB = config['contamination_check']['db']
    output:
        TSV = Path(config['working_dir']) / 'contamination_check' / 'kraken2' / 'tsv.io',
        TSV_report = Path(config['working_dir']) / 'contamination_check' / 'kraken2' / 'tsv-report.io',
        INFORMS = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'contamination_check' / 'kraken2'
    threads: 8
    priority: 1
    run:
        from camel.app.io.tooliodirectory import ToolIODirectory
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.kraken.kraken2 import Kraken2
        kraken2 = Kraken2(camel)
        kraken2.add_input_files(SnakePipelineUtils.extracts_fq_input(input.IO, key_pe='FASTQ_PE'))
        kraken2.add_input_files({'DB': [ToolIODirectory(input.DB)]})
        step = Step(rule, kraken2, camel, params.running_dir, config)
        kraken2.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(kraken2, output)

rule contamination_check_kraken_report_parser:
    """
    Parses the Kraken report and looks for contamination at the species level. 
    """
    input:
        TSV = rules.contamination_check_kraken2_run.output.TSV
    output:
        INFORMS = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'contamination_check' / 'kraken2',
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
        TSV = rules.contamination_check_kraken2_run.output.TSV
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
         INFORMS_kraken2 = rules.contamination_check_kraken2_run.output.INFORMS,
         TSV = rules.contamination_check_kraken2_run.output.TSV_report
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
            ('kraken2_expected_species', informs['expected'][0]),
            ('kraken2_expected_species_occurrence', informs['expected'][1]),
            ('kraken2_contaminants_warn', str(informs['contaminants_warn'])),
            ('kraken2_contaminants_fail', str(informs['contaminants_fail']))
        ]
        with open(output[0], 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')
