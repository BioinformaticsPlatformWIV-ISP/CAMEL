rule get_kraken_db:
    """
    Converts the Kraken database.
    """
    input:
        DB = config.get('db_kraken')
    output:
        os.path.join(__WORKING_DIR, 'contamination_check', 'db.io')
    run:
        print("Kraken DB: {}".format(input.DB))
        SnakemakeUtils.dump_object([ToolIODirectory(input.DB)], output[0])

rule kraken:
    """
    Kraken run to assign taxonomic labels to reads.
    """
    input:
        FASTQ_PE = os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-pe.io'),
        DB = os.path.join(__WORKING_DIR, 'contamination_check', 'db.io')
    output:
        TSV = os.path.join(__WORKING_DIR, 'contamination_check', 'tsv.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'contamination_check')
    run:
        from app.tools.kraken.kraken import Kraken
        kraken = Kraken(camel)
        SnakemakeUtils.add_pickle_inputs(kraken, input)
        step = SnakeStep(rule, kraken, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(kraken, output)

rule kraken_report:
    """
    Creates a report based on the Kraken output.
    """
    input:
        TSV = os.path.join(__WORKING_DIR, 'contamination_check', 'tsv.io'),
        DB = os.path.join(__WORKING_DIR, 'contamination_check', 'db.io')
    output:
        TSV = os.path.join(__WORKING_DIR, 'contamination_check', 'tsv-report.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'contamination_check')
    run:
        from app.tools.kraken.krakenreport import KrakenReport
        kraken_report = KrakenReport(camel)
        SnakemakeUtils.add_pickle_inputs(kraken_report, input)
        step = SnakeStep(rule, kraken_report, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(kraken_report, output)

rule kraken_report_parser:
    """
    Parses the Kraken report and looks for contamination at the species level.
    """
    input:
        TSV = os.path.join(__WORKING_DIR, 'contamination_check', 'tsv-report.io')
    output:
        INFORMS = os.path.join(__WORKING_DIR, 'contamination_check', 'informs-contamination.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'contamination_check')
    run:
        from app.tools.kraken.krakenreportparser import KrakenReportParser
        report_parser = KrakenReportParser(camel)
        # parameter updated in pipeline.step_tools_parameter table
        # report_parser.update_parameters(expected_species='Listeria monocytogenes')
        SnakemakeUtils.add_pickle_inputs(report_parser, input)
        step = SnakeStep(rule, report_parser, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(report_parser, output)

rule krona:
    """
    Creates an interactive pie chart displaying the Kraken output.
    """
    input:
        TSV = os.path.join(__WORKING_DIR, 'contamination_check', 'tsv.io')
    output:
        HTML = os.path.join(__WORKING_DIR, 'contamination_check', 'html-krona.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'contamination_check')
    run:
        from app.tools.krona.krona import Krona
        krona = Krona(camel)
        SnakemakeUtils.add_pickle_inputs(krona, input)
        step = SnakeStep(rule, krona, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(krona, output)

rule contamination_report:
    """
    Creates a report containing the results of the contamination check.
    """
    input:
        HTML_Krona = os.path.join(__WORKING_DIR, 'contamination_check', 'html-krona.io'),
        INFORMS_species = os.path.join(__WORKING_DIR, 'contamination_check', 'informs-contamination.io')
    output:
        VAL_HTML = os.path.join(__WORKING_DIR, 'contamination_check', 'html.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'contamination_check'),
        output_dir = config['output_dir']
    run:
        from app.tools.pipelines.quality_checks.htmlreportercontamination import HtmlReporterContamination
        reporter = HtmlReporterContamination(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = SnakeStep(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        reporter.tool_outputs['VAL_HTML'][0].value.copy_files(params.output_dir)
        SnakemakeUtils.dump_tool_outputs(reporter, output)
