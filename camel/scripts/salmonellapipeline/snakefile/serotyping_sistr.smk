import json
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import assembly
from camel.scripts.salmonellapipeline.snakefile import serotyping_sistr

camel = Camel.get_instance()

rule serotyping_sistr_run:
    """
    This rule executes SISTR to obtain serotpying results using cgMLST.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA
    output:
        JSON = Path(config['working_dir']) / 'serotyping_sistr' / 'sistr_output.io',
        INFORMS = Path(config['working_dir']) / 'serotyping_sistr' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping_sistr',
        db_path_sistr = config['serotyping']['sistr']['path']
    run:
        from camel.app.tools.pipelines.salmonella.sistr import Sistr

        sistr_tool = Sistr(camel)
        sistr_tool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path_sistr)))]})
        SnakemakeUtils.add_pickle_input(sistr_tool, 'FASTA', Path(input.FASTA))
        step = Step(str(rule), sistr_tool, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sistr_tool, output)

rule serotyping_sistr_dump_summary_info:
    """
    This rule creates a simple output summary for the SISTR serotyping analysis.
    """
    input:
        JSON_sistr = rules.serotyping_sistr_run.output.JSON,
        INFORMS_sistr = Path(config['working_dir']) / serotyping_sistr.OUTPUT_SEROTYPE_SISTR_INFORMS
    output:
        VAL_TSV_sistr = Path(config['working_dir']) / 'serotyping_sistr'/ 'summary_out_sistr.tsv'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping_sistr'
    run:
        import copy

        # parse SISTR output
        with SnakemakeUtils.load_object(Path(str(input.JSON_sistr)))[0].path.open('r') as handle:
            json_data = json.load(handle)[0]
        header_locus = ['Locus', 'serotype_or_group', '% Identity', 'HSP/Locus length', 'Contig', 'Position in contig']
        if json_data['qc_status'] == 'PASS':
            hits_dict_tsv = {
                'serotype_antigenic_formula':':'.join([
                    str(json_data['o_antigen']), str(json_data['h1']), str(json_data['h2'])]),
                'serotype_serogroup': json_data['serogroup'],
                'serotype_consensus': json_data['serovar'],
                'qc_status' : 'PASS'
            }
            hits_dict_json = copy.deepcopy(hits_dict_tsv)
            serotyping_sistr.sistr_output_parser(json_data['h1_flic_prediction'], 'fliC', 'h1', hits_dict_tsv)
            serotyping_sistr.sistr_output_parser(json_data['h2_fljb_prediction'], 'fljB', 'h2', hits_dict_tsv)
            serotyping_sistr.sistr_output_parser(json_data['serogroup_prediction']['wzx_prediction'], 'wzx', 'o', hits_dict_tsv)
            serotyping_sistr.sistr_output_parser(json_data['serogroup_prediction']['wzy_prediction'], 'wzy', 'o', hits_dict_tsv)
        else:
            hits_dict_tsv = {
                'serotype_antigenic_formula': '-',
                'serotype_serogroup': '-',
                'serotype_consensus': '-',
                'qc_status': 'FAIL'}
            for variable in ['hits_serotype_h1_fliC', 'hits_serotype_h2_fljB', 'hits_serotype_o_wzx', 'hits_serotype_o_wzy']:
                hits_dict_tsv[variable] = '-'

        informs_sistr = SnakemakeUtils.load_object(Path(str(input.INFORMS_sistr)))
        with Path(output.VAL_TSV_sistr).open('w') as handle:
            for k, v in hits_dict_tsv.items():
                line = f"sistr_{k}\t{v}\n"
                handle.write(line)
            handle.write(f"sistr_tool_version\t{informs_sistr['_name']}\n")
            handle.write(f"sistr_db_version\t{informs_sistr['last_update_date']}\n")

rule serotyping_sistr_report:
    """
    This rule creates a simple output report for the SISTR serotyping analysis.
    """
    input:
        JSON_SISTR = rules.serotyping_sistr_run.output.JSON,
        VAL_TSV = rules.serotyping_sistr_dump_summary_info.output.VAL_TSV_sistr,
        INFORMS_serotyping_sistr = rules.serotyping_sistr_run.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / 'serotyping_sistr' / 'html_sistr.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping_sistr' ,
        db_path_sistr = config['serotyping']['sistr']['path']
    run:
        from camel.app.tools.pipelines.salmonella.sistrreporter import SistrReporter

        reporter = SistrReporter(camel)
        reporter.add_input_files({'DIR_sistr': [ToolIODirectory(Path(str(params.db_path_sistr)))]})
        SnakemakeUtils.add_pickle_inputs(reporter, input, excluded_keys=['VAL_TSV'])
        reporter.add_input_files({'TSV_output': [ToolIOFile(Path(str(input.VAL_TSV)))]})
        step = Step(str(rule), reporter, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule serotyping_sistr_report_empty:
    """
    Creates an empty HTML report for the SISTR serotyping analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'serotyping_sistr' / 'html_sistr-empty.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping_sistr'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.pipelines.salmonella.sistrreporter import SistrReporter
        SnakePipelineUtils.create_empty_report_section(SistrReporter.TITLE, Path(output.VAL_HTML))