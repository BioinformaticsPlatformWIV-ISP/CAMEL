from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import abritamr, assembly


camel = Camel.get_instance()


rule abritamr_run:
    """
    This rule executes AbriTAMR run and gets the results.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA
    output:
        TXT_MATCHES = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_MATCHES,
        TXT_PARTIALS = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_PARTIALS,
        TSV_amrfinder = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_AMRFINDER,
        INFORMS = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_RUN_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'abritamr' ,
        db_path_amrf = config['abritamr']['amrfinderplus']['path'],
        species = config['abritamr']['species']
    run:
        from camel.app.tools.abritamr.abritamrrun import AbriTAMRRun

        abritamr_run_tool = AbriTAMRRun(camel)
        abritamr_run_tool.add_input_files({'DIR_AMRF': [ToolIODirectory(Path(str(params.db_path_amrf)))]})
        SnakemakeUtils.add_pickle_input(abritamr_run_tool, 'FASTA', Path(input.FASTA))
        step = Step(str(rule), abritamr_run_tool, camel, params.running_dir)
        abritamr_run_tool.update_parameters(species=params.species)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(abritamr_run_tool, output)

rule abritamr_generate_dummy_mdu_qc_file:
    """
    This rule creates a dummy QC file required by AbriTAMR report in order to make it work as expected.
    """
    output:
        TXT_MDU_QC = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_QC
    params:
        species = config['abritamr']['species'],
        running_dir = Path(config['working_dir']) / 'abritamr'
    run:
        species = params.species if not params.species == 'Salmonella' else 'Salmonella enterica'
        qc_file_txt = params.running_dir / 'qc_file.txt'
        with qc_file_txt.open('w') as handle:
            handle.write(f'ISOLATE,SPECIES_EXP,SPECIES_OBS,TEST_QC\n'
                         f'{params.running_dir},{species},{species},PASS\n')
        SnakemakeUtils.dump_object([ToolIOFile(qc_file_txt)], Path(output.TXT_MDU_QC))


rule abritamr_report_run:
    """
    This rule will run AbriTAMR report to generate an antibiogram in case the species is Salmonella, 
    else it will create a summary of the TXT_MATCHES and TXT_PARTIALS that will not actually be used in the 
    abritamr_report rule.
    """
    input:
        TXT_MDU_QC = rules.abritamr_generate_dummy_mdu_qc_file.output.TXT_MDU_QC,
        TXT_MATCHES = rules.abritamr_run.output.TXT_MATCHES,
        TXT_PARTIALS = rules.abritamr_run.output.TXT_PARTIALS,
        INFORMS_ABRITAMR_RUN = rules.abritamr_run.output.INFORMS
    output:
        REPORT_ABRITAMR = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_REPORT_REPORT,
        INFORMS = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_REPORT_REPORT_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'abritamr',
    run:
        from camel.app.tools.abritamr.abritamrreport import AbriTAMRReport

        abritamr_report_tool = AbriTAMRReport(camel)
        SnakemakeUtils.add_pickle_inputs(abritamr_report_tool, input)
        step = Step(str(rule), abritamr_report_tool, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(abritamr_report_tool, output)

rule abritamr_create_summary:
    """
    Creates the tabular summary output for the AbritAMR assay.
    """
    input:
        REPORT_ABRITAMR = rules.abritamr_report_run.output.REPORT_ABRITAMR,
        INFORMS_ABRITAMR_RUN = rules.abritamr_run.output.INFORMS,
        INFORMS_ABRITAMR_REPORT = rules.abritamr_report_run.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_SUMMARY
    run:
        import pandas as pd

        df_abritamr = pd.read_excel(
            SnakemakeUtils.load_object(Path(input.REPORT_ABRITAMR))[0].path, engine='openpyxl')
        df_abritamr.fillna('-', inplace=True)  # replace all missing values by dashes
        with Path(output.TSV).open('w') as handle:
            for column in df_abritamr.columns[2:]:
                key = f'abritamr_{column.replace(" - ","_")}'
                value = df_abritamr.iloc[0][column]
                handle.write(f"{key}\t{value}\n")

            informs_abritamr_run = SnakemakeUtils.load_object(Path(input.INFORMS_ABRITAMR_RUN))
            informs_abritamr_report = SnakemakeUtils.load_object(Path(input.INFORMS_ABRITAMR_REPORT))
            handle.write(f"abritamr_tool_version\tAbriTAMR {informs_abritamr_report['_version']}\n")
            handle.write(f"abritamr_db_version\t{informs_abritamr_run['last_update_date']}\n")

rule abritamr_report:
    """
    This rule creates a simple HTML output report for the AbriTAMR tool.
    """
    input:
        TSV = rules.abritamr_create_summary.output.TSV,
        REPORT_ABRITAMR = rules.abritamr_report_run.output.REPORT_ABRITAMR,
        INFORMS_ABRITAMR_RUN = rules.abritamr_run.output.INFORMS,
        TXT_MATCHES = rules.abritamr_run.output.TXT_MATCHES,
        TXT_PARTIALS = rules.abritamr_run.output.TXT_PARTIALS
    output:
        VAL_HTML = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'abritamr'
    run:
        from camel.app.tools.abritamr.abritamrreporter import AbriTAMRReporter

        reporter = AbriTAMRReporter(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input, excluded_keys=['TSV'])
        reporter.add_input_files({'TSV': [ToolIOFile(Path(input.TSV))]})
        step = Step(str(rule), reporter, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule abritamr_report_empty:
    """
    Creates an empty HTML report for the AbriTAMR analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section(AbriTAMRReporter.TITLE, Path(output.VAL_HTML))
