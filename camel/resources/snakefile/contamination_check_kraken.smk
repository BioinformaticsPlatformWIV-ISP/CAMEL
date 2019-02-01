"""
This Snakefile is used to perform a contamination check using KRAKEN. A report is generated with statistics and an
interactive KRONA visualization.
"""
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.contamination_check_kraken import OUTPUT_CONTAMINATION_SUMMARY, \
    OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY, OUTPUT_CONTAMINATION_CHECK_REPORT
from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_READS_PE
from camel.resources.snakefile.read_trimming_iontorrent import OUTPUT_TRIMMING_IT_READS

rule Contamination_check_get_db:
    """
    Converts the Kraken database to a pickle IO object.
    """
    output:
        os.path.join(config['working_dir'], 'contamination_check', 'db.io')
    params:
        db=config['contamination_check']['db']
    run:
        SnakemakeUtils.dump_object([ToolIODirectory(params.db)], output[0])

rule Contamination_check_kraken:
    """
    Assign taxonomic labels to reads using KRAKEN.
    """
    input:
        FASTQ_PE=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_PE) if config.get('read_type', 'illumina') == 'illumina' else [],
        FASTQ=os.path.join(config['working_dir'], OUTPUT_TRIMMING_IT_READS) if config.get('read_type', 'illumina') == 'iontorrent' else [],
        DB=os.path.join(config['working_dir'], 'contamination_check', 'db.io')
    output:
        TSV=os.path.join(config['working_dir'], 'contamination_check', 'kraken', 'tsv.io'),
        INFORMS=os.path.join(config['working_dir'], 'contamination_check', 'kraken', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'contamination_check', 'kraken')
    threads: 8
    priority: 1
    run:
        from camel.app.tools.kraken.kraken import Kraken
        kraken = Kraken(camel)
        SnakemakeUtils.add_pickle_inputs(kraken, input, ['DB'])
        for key in input.keys():
            if key == 'DB':
                continue
            if len(input[key]) > 0:
                kraken.add_input_files({key: SnakemakeUtils.load_object(input[key])})
        step = Step(rule, kraken, camel, params.running_dir, config)
        kraken.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(kraken, output)

rule Contamination_check_kraken_report:
    """
    Creteas a report based on the Kraken output.
    """
    input:
        TSV=os.path.join(config['working_dir'], 'contamination_check', 'kraken', 'tsv.io'),
        DB=os.path.join(config['working_dir'], 'contamination_check', 'db.io')
    output:
        TSV=os.path.join(config['working_dir'], 'contamination_check', 'kraken', 'tsv-report.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'contamination_check', 'kraken')
    run:
        from camel.app.tools.kraken.krakenreport import KrakenReport
        kraken_report = KrakenReport(camel)
        SnakemakeUtils.add_pickle_inputs(kraken_report, input)
        step = Step(rule, kraken_report, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(kraken_report, output)

rule Contamination_check_kraken_report_parser:
    """
    Parses the Kraken report and looks for contamination at the species level. 
    """
    input:
        TSV=os.path.join(config['working_dir'], 'contamination_check', 'kraken', 'tsv-report.io')
    output:
        INFORMS=os.path.join(config['working_dir'], 'contamination_check', 'kraken', 'informs-contamination.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'contamination_check', 'kraken'),
        expected_species=config['contamination_check']['expected_species']
    run:
        from camel.app.tools.kraken.krakenreportparser import KrakenReportParser
        report_parser = KrakenReportParser(camel)
        SnakemakeUtils.add_pickle_inputs(report_parser, input)
        step = Step(rule, report_parser, camel, params.running_dir, config)
        report_parser.update_parameters(expected_species=params.expected_species)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(report_parser, output)

rule Contamination_check_krona:
    """
    Creates an interactive pie chart displaying the Kraken output.
    """
    input:
        TSV=os.path.join(config['working_dir'], 'contamination_check', 'kraken', 'tsv.io')
    output:
        HTML=os.path.join(config['working_dir'], 'contamination_check', 'krona', 'html.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'contamination_check', 'krona')
    run:
        from camel.app.tools.krona.krona import Krona
        krona = Krona(camel)
        SnakemakeUtils.add_pickle_inputs(krona, input)
        step = Step(rule, krona, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(krona, output)

rule Contamination_check_report:
    """
    Creates a report containing the results of the contamination check.
    """
    input:
        HTML_Krona=os.path.join(config['working_dir'], 'contamination_check', 'krona', 'html.io'),
        INFORMS_species=os.path.join(config['working_dir'], 'contamination_check', 'kraken', 'informs-contamination.io'),
        INFORMS_kraken=os.path.join(config['working_dir'], 'contamination_check', 'kraken', 'informs.io'),
        TSV=os.path.join(config['working_dir'], 'contamination_check', 'kraken', 'tsv-report.io'),
        DB=os.path.join(config['working_dir'], 'contamination_check', 'db.io')
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_CONTAMINATION_CHECK_REPORT)
    params:
        running_dir=os.path.join(config['working_dir'], 'contamination_check', 'report')
    run:
        from camel.app.tools.pipelines.quality_checks.htmlreportercontamination import HtmlReporterContamination
        reporter = HtmlReporterContamination(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule Contamination_check_report_empty:
    """
    Generates an empty contamination check report.
    """
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY)
    params:
        running_dir=os.path.join(config['working_dir'], 'contamination_check', 'report')
    run:
        from camel.app.tools.pipelines.quality_checks.htmlreportercontamination import HtmlReporterContamination
        section = HtmlReporterContamination.generate_empty_section()
        SnakemakeUtils.dump_object([ToolIOValue(section)], output[0])

rule Contamination_check_dump_summary_info:
    """
    Dumps the summary information for the contamination check in tabular format.
    """
    input:
        INFORMS_species=os.path.join(config['working_dir'], 'contamination_check', 'kraken', 'informs-contamination.io')
    output:
        os.path.join(config['working_dir'], OUTPUT_CONTAMINATION_SUMMARY)
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
