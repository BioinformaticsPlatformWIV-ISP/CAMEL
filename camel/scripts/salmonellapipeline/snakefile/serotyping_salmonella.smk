import json
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import assembly
from camel.scripts.salmonellapipeline.snakefile import serotyping_salmonella

camel = Camel.get_instance()

rule serotyping_sistr_run:
    """
    This rule executes SISTR to obtain serotpying results using cgMLST.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA
    output:
        JSON = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'serotyping_sistr' / 'sistr_output.io',
        INFORMS = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'serotyping_sistr' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / wildcards.input_format / 'serotyping_sistr',
        db_path_sistr = config['serotyping']['sistr']['path']
    run:
        from camel.app.tools.pipelines.salmonella.sistr import Sistr

        sistr_tool = Sistr(camel)
        sistr_tool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path_sistr)))]})
        SnakemakeUtils.add_pickle_input(sistr_tool, 'FASTA', Path(input.FASTA))
        step = Step(str(rule), sistr_tool, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sistr_tool, output)

rule serotyping_seqsero2_run:
    """
    This rule executes SeqSero2 in the appropriate modes depending on if input is fasta or not.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA,
        IO = Path(config['working_dir']) / 'fq_dict.io' if 'fasta' not in config['input'] else []
    output:
        TXT = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'serotyping_seqsero2_{mode}' / 'SeqSero2_result.io',
        INFORMS = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'serotyping_seqsero2_{mode}' / 'informs.io'
    params:
        mode = lambda wildcards: wildcards.mode,
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / wildcards.input_format / f'serotyping_seqsero2_{wildcards.mode}',
        db_path_seqsero2 = config['serotyping']['seqsero2']['path'],
        read_key = lambda wildcards: wildcards.input_format
    run:
        from camel.app.tools.pipelines.salmonella.seqsero2 import SeqSero2

        seqsero_tool = SeqSero2(camel)
        seqsero_tool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path_seqsero2)))]})
        SnakemakeUtils.add_pickle_input(seqsero_tool, 'FASTA', Path(input.FASTA))
        if params.read_key == 'fastq_pe':
            seqsero_tool.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        if params.read_key == 'fastq_se':
            if config['input_type'] == 'ont':
                seqsero_tool.add_input_files(SnakePipelineUtils.extracts_fq_input(
                    Path(input.IO),key_se='FASTQ_ONT',read_type='SE'))
            else:
                seqsero_tool.add_input_files(SnakePipelineUtils.extracts_fq_input(
                    Path(input.IO), key_se='FASTQ', read_type='SE'))
        seqsero_tool.update_parameters(mode=str(params.mode))
        step = Step(str(rule), seqsero_tool, camel, Path(str(params.running_dir)))
        step.run_step()
        if config['input_type'] == 'hybrid':
            seqsero_tool.informs['_tag'] = f"{params.mode} - {'Illumina' if params.read_key == 'fastq_pe' else 'ONT'}"
        SnakemakeUtils.dump_tool_outputs(seqsero_tool, output)

rule serotyping_sistr_dump_summary_info:
    """
    This rule creates a simple output summary for the SISTR serotyping analysis.
    """
    input:
        JSON_sistr = rules.serotyping_sistr_run.output.JSON,
        INFORMS_sistr = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEROTYPE_SISTR_INFORMS
    output:
        VAL_TSV_sistr = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'summary_out_sistr.tsv'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / wildcards.input_format
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
            serotyping_salmonella.sistr_output_parser(json_data['h1_flic_prediction'], 'fliC', 'h1', hits_dict_tsv)
            serotyping_salmonella.sistr_output_parser(json_data['h2_fljb_prediction'], 'fljB', 'h2', hits_dict_tsv)
            serotyping_salmonella.sistr_output_parser(json_data['serogroup_prediction']['wzx_prediction'], 'wzx', 'o', hits_dict_tsv)
            serotyping_salmonella.sistr_output_parser(json_data['serogroup_prediction']['wzy_prediction'], 'wzy', 'o', hits_dict_tsv)
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

rule serotyping_seqsero2_dump_summary_info:
    """
    This rule creates a simple output summary for the SeqSero2 serotyping analysis.
    """
    input:
        TXT_seqsero2_kmer = lambda wildcards: str(rules.serotyping_seqsero2_run.output.TXT).format(mode='Kmer', input_format=wildcards.input_format),
        INFORMS_seqsero2_kmer = lambda wildcards: str(rules.serotyping_seqsero2_run.output.INFORMS).format(mode='Kmer', input_format=wildcards.input_format),
        TXT_seqsero2_allele = lambda wildcards: str(rules.serotyping_seqsero2_run.output.TXT).format(mode='Allele', input_format=wildcards.input_format) if 'fasta' not in config['input'] else [],
        INFORMS_seqsero2_allele = lambda wildcards: str(rules.serotyping_seqsero2_run.output.INFORMS).format(mode='Allele', input_format=wildcards.input_format) if 'fasta' not in config['input'] else [],
        TXT_seqsero2_kmerread = lambda wildcards: str(rules.serotyping_seqsero2_run.output.TXT).format(mode='Kmerread', input_format=wildcards.input_format) if 'fasta' not in config['input'] else [],
        INFORMS_seqsero2_kmerread = lambda wildcards: str(rules.serotyping_seqsero2_run.output.INFORMS).format(mode='Kmerread', input_format=wildcards.input_format) if 'fasta' not in config['input'] else []
    output:
        VAL_TSV_seqsero2 = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'summary_out_seqsero2.tsv'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / wildcards.input_format
    run:
        # parse obligate SeqSero2 output
        informs_seqsero2_kmer = SnakemakeUtils.load_object(Path(str(input.INFORMS_seqsero2_kmer)))
        tsv_results = serotyping_salmonella.seqsero2_output_parser(SnakemakeUtils.load_object(Path(str(input.TXT_seqsero2_kmer)))[0].path, 'seqsero2_kmer')
        with Path(output.VAL_TSV_seqsero2).open('w') as handle:
            handle.writelines(item + '\n' for item in tsv_results)

        # parse facultative SeqSero2 output
        if 'fasta' not in config['input']:
            for args_tuple in [(SnakemakeUtils.load_object(Path(str(input.TXT_seqsero2_allele)))[0].path, 'seqsero2_allele'),
                               (SnakemakeUtils.load_object(Path(str(input.TXT_seqsero2_kmerread)))[0].path, 'seqsero2_kmerread')
                                ]:
                tsv_results = serotyping_salmonella.seqsero2_output_parser(*args_tuple)
                with Path(output.VAL_TSV_seqsero2).open('a') as handle:
                    for item in tsv_results:
                        handle.write(item + '\n')

        with Path(output.VAL_TSV_seqsero2).open('a') as handle:
            handle.write(f"seqsero2_tool_version\t{informs_seqsero2_kmer['_name']}\n")
            handle.write(f"seqsero2_db_version\t{informs_seqsero2_kmer['last_update_date']}\n")

rule serotyping_sistr_report:
    """
    This rule creates a simple output report for the SISTR serotyping analysis.
    """
    input:
        JSON_SISTR = rules.serotyping_sistr_run.output.JSON,
        VAL_TSV = rules.serotyping_sistr_dump_summary_info.output.VAL_TSV_sistr,
        INFORMS_serotyping_sistr = rules.serotyping_sistr_run.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'html_sistr.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / wildcards.input_format,
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

rule serotyping_seqsero2_report:
    """
    This rule creates a simple output report for the SeqSero2 serotyping analysis.
    """
    input:
        TXT_seqsero2_kmer = lambda wildcards: str(rules.serotyping_seqsero2_run.output.TXT).format(mode='Kmer', input_format=wildcards.input_format),
        TXT_seqsero2_allele = lambda wildcards: str(rules.serotyping_seqsero2_run.output.TXT).format(mode='Allele', input_format=wildcards.input_format) if 'fasta' not in config['input'] else [],
        TXT_seqsero2_kmerread = lambda wildcards: str(rules.serotyping_seqsero2_run.output.TXT).format(mode='Kmerread', input_format=wildcards.input_format) if 'fasta' not in config['input'] else [],
        VAL_TSV = rules.serotyping_seqsero2_dump_summary_info.output.VAL_TSV_seqsero2,
        INFORMS_serotyping_seqsero2 = lambda wildcards: str(rules.serotyping_seqsero2_run.output.INFORMS).format(mode='Kmer', input_format=wildcards.input_format)
    output:
        VAL_HTML = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'html_seqsero2.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / wildcards.input_format,
        db_path_seqsero2 = config['serotyping']['seqsero2']['path']
    run:
        from camel.app.tools.pipelines.salmonella.seqsero2reporter import SeqSero2Reporter

        reporter = SeqSero2Reporter(camel)
        reporter.add_input_files({'DIR_seqsero2': [ToolIODirectory(Path(str(params.db_path_seqsero2)))]})
        SnakemakeUtils.add_pickle_inputs(reporter, input, excluded_keys=['VAL_TSV', 'TXT_seqsero2_allele', 'TXT_seqsero2_kmerread'])
        reporter.add_input_files({'TSV_output': [ToolIOFile(Path(input.VAL_TSV))]})
        if input.TXT_seqsero2_allele:
            SnakemakeUtils.add_pickle_input(reporter, 'TXT_seqsero2_allele', Path(str(input.TXT_seqsero2_allele)))
        if input.TXT_seqsero2_kmerread:
            SnakemakeUtils.add_pickle_input(reporter, 'TXT_seqsero2_kmerread', Path(str(input.TXT_seqsero2_kmerread)))
        step = Step(str(rule), reporter, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule serotyping_sistr_report_empty:
    """
    Creates an empty HTML report for the SISTR serotyping analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'html_sistr-empty.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / wildcards.input_format
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.pipelines.salmonella.sistrreporter import SistrReporter
        SnakePipelineUtils.create_empty_report_section(SistrReporter.TITLE, Path(output.VAL_HTML))

rule serotyping_seqsero2_report_empty:
    """
    Creates an empty HTML report for the SeqSero2 serotyping analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'html_seqsero2-empty.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / wildcards.input_format
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.pipelines.salmonella.seqsero2reporter import SeqSero2Reporter
        SnakePipelineUtils.create_empty_report_section(SeqSero2Reporter.TITLE, Path(output.VAL_HTML))
