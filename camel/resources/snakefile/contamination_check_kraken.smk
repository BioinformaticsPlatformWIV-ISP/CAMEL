from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import contamination_check_kraken


rule contamination_check_kraken2_run:
    """
    Assigns taxonomic labels to reads using KRAKEN2.
    """
    input:
        IO = str(contamination_check_kraken.get_input(config)),
        DB = config['contamination_check']['db']
    output:
        TSV = Path(config['working_dir']) / 'contamination_check' / '{input_format}' / 'kraken2' / 'tsv.io',
        TSV_report = Path(config['working_dir']) / 'contamination_check' / '{input_format}' / 'kraken2' / 'tsv-report.io',
        INFORMS = Path(config['working_dir']) / 'contamination_check' / '{input_format}' / 'kraken2' / 'informs.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'contamination_check' / wildcards.input_format / 'kraken2',
        input_format = lambda wildcards: wildcards.input_format
    threads: 8
    priority: 1
    run:
        from camel.app.io.tooliodirectory import ToolIODirectory
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.kraken.kraken2 import Kraken2

        kraken2 = Kraken2(Camel.get_instance())
        if params.input_format == 'fastq_pe':
            kraken2.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        elif params.input_format == 'fastq_se':
            kraken2.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_se='FASTQ', read_type='SE'))
        else:
            SnakemakeUtils.add_pickle_input(kraken2, 'FASTA', Path(input.IO))
        kraken2.add_input_files({'DB': [ToolIODirectory(Path(input.DB))]})
        step = Step(str(rule), kraken2, Camel.get_instance(), Path(str(params.dir_)))
        kraken2.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(kraken2, output)

rule contamination_check_kraken_report_parser:
    """
    Parses the Kraken report and looks for contamination at the species level. 
    """
    input:
        TSV = rules.contamination_check_kraken2_run.output.TSV_report
    output:
        INFORMS = Path(config['working_dir']) / 'contamination_check' / '{input_format}' / 'kraken2' / 'informs-contamination.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'contamination_check' / wildcards.input_format / 'kraken2',
        expected_species = config['contamination_check']['expected_species'],
        allowed_species = config['contamination_check'].get('allowed_species', None),
        level_of_depth = config['contamination_check'].get('level_of_depth', 'S')
    run:
        from camel.app.tools.kraken.krakenreportparser import KrakenReportParser
        report_parser = KrakenReportParser(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(report_parser, input)
        step = Step(str(rule), report_parser, Camel.get_instance(), Path(str(params.dir_)))
        report_parser.update_parameters(expected_species=params.expected_species)
        if params.level_of_depth is not None:
            report_parser.update_parameters(level_of_depth=params.level_of_depth)
        if params.allowed_species is not None:
            report_parser.update_parameters(allowed_species=','.join(params.allowed_species))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(report_parser, output)

rule contamination_check_krona:
    """
    Creates an interactive pie chart displaying the Kraken output.
    """
    input:
        TSV = rules.contamination_check_kraken2_run.output.TSV
    output:
        HTML = Path(config['working_dir']) / 'contamination_check' / '{input_format}' / 'krona' / 'html.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'contamination_check' / wildcards.input_format / 'krona'
    run:
        from camel.app.tools.krona.krona import Krona
        krona = Krona(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(krona, input)
        step = Step(str(rule), krona, Camel.get_instance(), Path(str(params.dir_)))
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
        VAL_HTML = Path(config['working_dir']) / 'contamination_check' / '{input_format}' / 'report' / 'html.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'contamination_check' / wildcards.input_format / 'report',
        input_format = lambda wildcards: wildcards.input_format,
        input_type = config['input_type']
    run:
        from camel.app.tools.pipelines.quality_checks.htmlreportercontamination import HtmlReporterContamination
        reporter = HtmlReporterContamination(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)

        # Update report title + output files for hybrid data
        if params.input_type == 'hybrid':
            reporter.update_parameters(
                suffix='illumina' if params.input_format == 'fastq_pe' else 'ont',
                suffix_title='Illumina' if params.input_format == 'fastq_pe' else 'ONT',
            )

        # Run the tool
        step = Step(str(rule), reporter, Camel.get_instance(), Path(str(params.dir_)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule contamination_check_report_empty:
    """
    Generates an empty contamination check report.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'contamination_check' / '{input_format}' / 'report' / 'html-empty.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'contamination_check' / wildcards.input_format / 'report',
        input_format = lambda wildcards: wildcards.input_format,
        input_type = config['input_type']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.quality_checks.htmlreportercontamination import HtmlReporterContamination

        # Suffix for hybrid data
        if params.input_type == 'hybrid':
            suffix_title = 'Illumina' if params.input_format == 'fastq_pe' else 'ONT'
        else:
            suffix_title = None
        section = HtmlReporterContamination.generate_empty_section(suffix_title)
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.VAL_HTML))

rule contamination_check_dump_summary_info:
    """
    Dumps the summary information for the contamination check in tabular format.
    """
    input:
        INFORMS_kraken = rules.contamination_check_kraken2_run.output.INFORMS,
        INFORMS_species = rules.contamination_check_kraken_report_parser.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / 'contamination_check' / '{input_format}' / 'summary' / 'summary_out.tsv'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'contamination_check' / wildcards.input_format / 'summary',
        input_type = config['input_type'],
        input_format = lambda wildcards: wildcards.input_format,
        allowed_species = config['contamination_check'].get('allowed_species', None)
    run:
        informs_kraken2 = SnakemakeUtils.load_object(Path(input.INFORMS_kraken))
        informs = SnakemakeUtils.load_object(Path(input.INFORMS_species))

        suffix = f'_{params.input_format}' if params.input_type == 'hybrid' else ''
        summary_data = [
            (f'kraken2{suffix}_expected_taxon', informs['expected'][0]),
            (f'kraken2{suffix}_expected_taxon_occurrence', informs['expected'][1]),
            (f'kraken2{suffix}_contaminants_warn', str(informs['contaminants_warn'])),
            (f'kraken2{suffix}_contaminants_fail', str(informs['contaminants_fail'])),
        ]
        if params.allowed_species is not None:
            summary_data.append((f'kraken2{suffix}_allowed', str(informs['allowed'])))
        summary_data.extend([
            (f'kraken2{suffix}_tool_version', informs_kraken2['_name']),
            (f'kraken2{suffix}_db', informs_kraken2['database']['name']),
            (f'kraken2{suffix}_last_update', informs_kraken2['database']['last_update'])
        ])
        with open(output.TSV, 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')
