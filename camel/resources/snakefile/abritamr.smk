from pathlib import Path

import numpy as np

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.tools.abritamr.abritamrreport import AbriTAMRReport
from camel.app.tools.abritamr.abritamrreporter import AbriTAMRReporter
from camel.app.tools.abritamr.abritamrrun import AbriTAMRRun
from camel.resources.snakefile import abritamr, assembly


camel = Camel.get_instance()


rule abritamr_run:
    """
    This rule executes abritamr run and gets the results
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA
    output:
        TXT_MATCHES = Path(config['working_dir']) /  abritamr.OUTPUT_MATCHES_ABRITAMR,
        TXT_PARTIALS = Path(config['working_dir']) / abritamr.OUTPUT_PARTIALS_ABRITAMR,
        AMRFINDER_OUT = Path(config['working_dir']) / abritamr.OUTPUT_AMRFINDER_ABRITAMR,
        INFORMS = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_RUN_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'abritamr' ,
        db_path_amrf = config['abritamr']['amrfinderplus']['path'],
        species = config['abritamr']['species'],
        jobs = 4
    run:
        abritamrruntool = AbriTAMRRun(camel)
        abritamrruntool.add_input_files({'DIR_AMRF': [ToolIODirectory(Path(str(params.db_path_amrf)))]})
        SnakemakeUtils.add_pickle_input(abritamrruntool,'FASTA',Path(input.FASTA))
        step = Step(str(rule),abritamrruntool,camel,params.running_dir,config)
        abritamrruntool.update_parameters(jobs=params.jobs,species=params.species)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(abritamrruntool,output)

rule generate_dummy_mdu_qc_file:
    """
    This rule creates a qc file required by abritamr report in order to make it work as expected
    """
    output:
        TXT_MDU_QC = Path(config['working_dir']) /  abritamr.OUTPUT_QC_ABRITAMR
    params:
        species = config['abritamr']['species'],
        running_dir = Path(config['working_dir']) / 'abritamr'
    run:
        species = params.species if not params.species == 'Salmonella' else 'Salmonella enterica'
        with Path(output.TXT_MDU_QC).open('w') as handle:
            handle.write(f'ISOLATE,SPECIES_EXP,SPECIES_OBS,TEST_QC\n{params.running_dir},{species},{species},PASS\n')

rule abritamr_report_run:
    """
    This rule will run abritamr report to generate the final files
    """
    input:
        TXT_MDU_QC = rules.generate_dummy_mdu_qc_file.output.TXT_MDU_QC,
        TXT_MATCHES = rules.abritamr_run.output.TXT_MATCHES,
        TXT_PARTIALS = rules.abritamr_run.output.TXT_PARTIALS
    output:
        REPORT_ABRITAMR = Path(config['working_dir']) / abritamr.OUTPUT_REPORT_ABRITAMR,
        INFORMS = Path(config['working_dir']) / abritamr.OUTPUT_REPORT_ABRITAMR_INFORMS
    params:
        running_dir=Path(config['working_dir']) / 'abritamr',
        species = config['abritamr']['species']
    run:
        abritamrreporttool = AbriTAMRReport(camel)
        SnakemakeUtils.add_pickle_input(abritamrreporttool, 'TXT_MATCHES', Path(input.TXT_MATCHES))
        SnakemakeUtils.add_pickle_input(abritamrreporttool, 'TXT_PARTIALS', Path(input.TXT_PARTIALS))
        abritamrreporttool.add_input_files({'TXT_MDU_QC': [ToolIOFile(Path(input.TXT_MDU_QC))],
                                            'VAL_SPECIES': [ToolIOValue(params.species)]})
        step = Step(str(rule), abritamrreporttool, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(abritamrreporttool, output)

rule create_output_summary_abritamr:
    """
    This rule creates the output summary to generate the summary tsv and the json for the tool results
    """
    input:
        REPORT_ABRITAMR = rules.abritamr_report_run.output.REPORT_ABRITAMR,
        INFORMS_ABRITAMR_RUN = rules.abritamr_run.output.INFORMS,
        INFORMS_ABRITAMR_REPORT = rules.abritamr_report_run.output.INFORMS
    output:
        VAL_TSV = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_SUMMARY,
        JSON = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_SUMMARY_JSON
    params:
        species = config['abritamr']['species']
    run:
        import pandas, json

        json_dict = {}
        informs_abritamr_run = SnakemakeUtils.load_object(Path(input.INFORMS_ABRITAMR_RUN))
        informs_abritamr_report = SnakemakeUtils.load_object(Path(input.INFORMS_ABRITAMR_REPORT))
        df_abritamr = pandas.read_excel(
            SnakemakeUtils.load_object(Path(input.REPORT_ABRITAMR))[0].path,engine='openpyxl')
        #replace all nan by dashes
        df_abritamr.replace(np.nan,'-',inplace=True)
        with Path(output.VAL_TSV).open('w') as file_tsv:
            for column in range(2,len(df_abritamr.columns)):
                key = f'abritamr_{df_abritamr.columns[column].replace(" - ","_")}'
                value = df_abritamr.iloc[0, column]
                file_tsv.write(f"{key}\t{value}\n")
                json_dict[key] = value
            file_tsv.write(f"abritamr_tool_version\tAbriTAMR {informs_abritamr_report['_version']}\n")
            file_tsv.write(f"abritamr_db_version\t{informs_abritamr_run['last_update_date']}\n")
        meta_json_dict = {'abritamr': {'results': json_dict,
                                       'informs_tools': {
                                           informs_abritamr_report.get('_tool',informs_abritamr_report['_name']): {
                                               '_name': informs_abritamr_report['_name'],
                                               '_version': informs_abritamr_report['_version'],
                                               '_command': informs_abritamr_report['_command']},
                                           informs_abritamr_run.get('_tool',informs_abritamr_run['_name']): {
                                               '_name': informs_abritamr_run['_name'],
                                               '_version': informs_abritamr_run['_version'],
                                               '_command': informs_abritamr_run['_command']}
                                       },
                                       'informs_dbs': {'last_updated': informs_abritamr_run['last_update_date'],
                                                       'name': informs_abritamr_run['key'],
                                                       'title': informs_abritamr_run['key']}
                                       }
                          }
        with Path(output.JSON).open('w') as handle:
            handle.write(json.dumps(meta_json_dict))

rule create_output_report_abritamr:
    """
    This rule creates a simple html output report for the AbriTAMR tool.
    """
    input:
        VAL_TSV = rules.create_output_summary_abritamr.output.VAL_TSV,
        INFORMS_ABRITAMR_RUN = rules.abritamr_run.output.INFORMS,
        TXT_MATCHES = rules.abritamr_run.output.TXT_MATCHES,
        TXT_PARTIALS = rules.abritamr_run.output.TXT_PARTIALS
    output:
        VAL_HTML = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_REPORT
    params:
        species = config['abritamr']['species'],
        running_dir = Path(config['working_dir']) / 'abritamr'
    run:
        reportertool = AbriTAMRReporter(camel)
        SnakemakeUtils.add_pickle_inputs(reportertool, input, excluded_keys=['VAL_TSV'])
        reportertool.add_input_files({'TSV_output': [ToolIOFile(Path(input.VAL_TSV))], 'VAL_SPECIES': [ToolIOValue(params.species)]})
        step = Step(str(rule), reportertool, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reportertool, output)

rule abritamr_report_empty:
    """
    Creates an empty HTML report for the AbriTAMR analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / abritamr.OUTPUT_ABRITAMR_REPORT_EMPTY
    params:
        running_dir = Path(config['working_dir']) / 'abritamr'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section(AbriTAMRReporter.TITLE, Path(output.VAL_HTML))

