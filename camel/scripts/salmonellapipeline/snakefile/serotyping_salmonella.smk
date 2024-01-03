import json
from pathlib import Path
from typing import Any, Dict

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.pipelines.salmonella.seqsero2 import SeqSero2
from camel.app.tools.pipelines.salmonella.serovarsalmonellareporter import SerovarSalmonellaReporter
from camel.app.tools.pipelines.salmonella.sistr import Sistr
from camel.resources.snakefile import assembly_spades
from camel.scripts.salmonellapipeline.snakefile import serotyping_salmonella

camel = Camel.get_instance()

rule serotyping_sistr:
    """
    This rule executes SISTR to obtain serotpying results using cgMLST.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA
    output:
        JSON = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SISTR_SEROTYPE,
        INFORMS = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEROTYPE_SISTR_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'serotyping' / 'serotyping_sistr',
        db_path_sistr = config['serotyping']['sistr']['path']
    run:
        sistrtool = Sistr(camel)
        sistrtool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path_sistr)))]})
        SnakemakeUtils.add_pickle_input(sistrtool, 'FASTA', Path(input.FASTA))
        step = Step(str(rule), sistrtool, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sistrtool, output)

rule serotyping_seqsero2_wildcards:
    """
    This rule executes SeqSero2 in kmer mode.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA,  # todo change io to fasta input from anywhere (e.g. human read scrubbing) using link fasta in main.smk
        IO = Path(config['working_dir']) / 'fq_dict.io' if 'fasta' not in config['input'] else []  # todo change io to fastq input from anywhere (e.g. human read scrubbing) using link fastq in main.smk
    output:
        TXT = Path(config['working_dir']) / 'serotyping' / '_'.join(['serotyping_seqsero2', '{mode}']) / 'SeqSero2_result.io',
        INFORMS = Path(config['working_dir']) / 'serotyping' / '_'.join(['serotyping_seqsero2', '{mode}']) / 'informs.io'
    params:
        seqsero2_mode = lambda wildcards: wildcards.mode,
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / f'serotyping_seqsero2_{wildcards.mode}',
        db_path_seqsero2 = config['serotyping']['seqsero2']['path'],
        read_type= 'SE' if config.get('read_type') == 'iontorrent' else 'PE'
    run:
        seqserotool = SeqSero2(camel)
        seqserotool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path_seqsero2)))], 'MODE': [ToolIOValue(str(params.seqsero2_mode))]})
        SnakemakeUtils.add_pickle_input(seqserotool, 'FASTA', Path(input.FASTA))
        if params.read_type == 'PE':
            seqserotool.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        else:
            seqserotool.add_input_files(SnakePipelineUtils.extracts_fq_input(
                Path(input.IO), key_se='FASTQ', read_type=params.read_type))
        step = Step(str(rule), seqserotool, camel, Path(str(params.running_dir)), config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqserotool, output)

rule serotyping_dump_summary_info:
    """
    This rule creates a simple output report for SISTR and SeqSero2.
    """
    input:
        JSON_sistr = rules.serotyping_sistr.output.JSON,
        INFORMS_sistr = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEROTYPE_SISTR_INFORMS,
        TXT_seqsero2_kmer = expand(rules.serotyping_seqsero2_wildcards.output.TXT, mode='Kmer'),
        INFORMS_seqsero2_kmer = expand(rules.serotyping_seqsero2_wildcards.output.INFORMS, mode='Kmer'),
        TXT_seqsero2_allele = expand(rules.serotyping_seqsero2_wildcards.output.TXT, mode='Allele') if 'fasta' not in config['input'] else [],
        INFORMS_seqsero2_allele = expand(rules.serotyping_seqsero2_wildcards.output.INFORMS, mode='Allele') if 'fasta' not in config['input'] else [],
        TXT_seqsero2_kmerread = expand(rules.serotyping_seqsero2_wildcards.output.TXT, mode='Kmerread') if 'fasta' not in config['input'] else [],
        INFORMS_seqsero2_kmerread = expand(rules.serotyping_seqsero2_wildcards.output.INFORMS, mode='Kmerread') if 'fasta' not in config['input'] else []
    output:
        VAL_TSV = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEROTYPE_SUMMARY,
        JSON = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEROTYPE_SUMMARY_JSON
    params:
        running_dir = Path(config['working_dir']) / 'serotyping'
    threads: 8
    run:
        import copy
        meta_json_dict = {}

        # parse obligate Sistr output
        with SnakemakeUtils.load_object(Path(input.JSON_sistr))[0].path.open('r') as handle:
            json_data = json.load(handle)[0]
        header_locus = ['Locus', 'serotype_or_group', '% Identity', 'HSP/Locus length', 'Contig', 'Position in contig']
        if json_data['qc_status'] == 'PASS':
            hits_dict_tsv = {'serotype_antigenic_formula':':'.join([str(json_data['o_antigen']),
                                                                    str(json_data['h1']),
                                                                    str(json_data['h2'])
                                                                    ]),
                             'serotype_serogroup': json_data['serogroup'],
                             'serotype_concensus': json_data['serovar'],
                             'qc_status' : 'PASS'
                             }
            hits_dict_json = copy.deepcopy(hits_dict_tsv)
            serotyping_salmonella.sistr_output_parser(json_data['h1_flic_prediction'], 'fliC', 'h1', hits_dict_tsv, hits_dict_json, header_locus)
            serotyping_salmonella.sistr_output_parser(json_data['h2_fljb_prediction'], 'fljB', 'h2', hits_dict_tsv, hits_dict_json, header_locus)
            serotyping_salmonella.sistr_output_parser(json_data['serogroup_prediction']['wzx_prediction'], 'wzx', 'o', hits_dict_tsv, hits_dict_json, header_locus)
            serotyping_salmonella.sistr_output_parser(json_data['serogroup_prediction']['wzy_prediction'], 'wzy', 'o', hits_dict_tsv, hits_dict_json, header_locus)
        else:
            hits_dict_tsv = {'serotype_antigenic_formula': '-',
                            'serotype_serogroup': '-',
                            'serotype_concensus': '-',
                             'qc_status': 'FAIL'
                             }
            hits_dict_json: Dict[str, Any] = copy.deepcopy(hits_dict_tsv)
            for variable in ['hits_serotype_h1_fliC', 'hits_serotype_h2_fljB', 'hits_serotype_o_wzx', 'hits_serotype_o_wzy']:
                hits_dict_tsv[variable] = '-'
                hits_dict_json[variable] = {item: '-' for item in header_locus}

            with output.VAL_TSV.open('w') as handle:
                handle.writelines('\t'.join([f"sistr_{k}", v]) + '\n' for k, v in hits_dict_tsv.items())
            informs_sistr = SnakemakeUtils.load_object(Path(input.INFORMS_sistr))
            meta_json_dict.update({'sistr' : {**hits_dict_json, 'informs_tools' : { informs_sistr['_tool']: {'_name': informs_sistr['_name'], '_version': informs_sistr['_version'], '_command': informs_sistr['_command']}}, 'informs_dbs' : {'last_updated': informs_sistr['last_update_date'], 'name': informs_sistr['key'], 'title': informs_sistr['key']} }})

        # parse obligate seqsero2 output
        inter_json_dict, tsv_results = serotyping_salmonella.seqsero2_output_parser(SnakemakeUtils.load_object(Path(input.TXT_seqsero2_kmer))[0].path, 'seqsero2_kmer', SnakemakeUtils.load_object(Path(input.INFORMS_seqsero2_kmer)))
        meta_json_dict.update(inter_json_dict)
        with output.VAL_TSV.open('a') as handle:
            handle.writelines(item + '\n' for item in tsv_results)

        # parse facultative seqsero2 output
        if 'fasta' not in config['input']:
            for args_tuple in [(SnakemakeUtils.load_object(Path(input.TXT_seqsero2_allele))[0].path, 'seqsero2_allele', SnakemakeUtils.load_object(Path(input.INFORMS_seqsero2_allele))),
                               (SnakemakeUtils.load_object(Path(input.TXT_seqsero2_kmerread))[0].path, 'seqsero2_kmerread', SnakemakeUtils.load_object(Path(input.INFORMS_seqsero2_kmerread)))
                                ]:
                inter_json_dict, tsv_results = serotyping_salmonella.seqsero2_output_parser(args_tuple[0], args_tuple[1], args_tuple[2])
                meta_json_dict.update(inter_json_dict)
                with output.VAL_TSV.open('a') as handle:
                    handle.writelines(item + '\n' for item in tsv_results)

        with output.JSON.open('w') as handle:
            handle.write(json.dumps(meta_json_dict))

rule create_output_report_serotyping:
    """
    This rule creates a simple output report, combining both serotyping tools
    """
    input:
        TXT_seqsero2_kmer = expand(rules.serotyping_seqsero2_wildcards.output.TXT, mode='Kmer'),
        TXT_seqsero2_allele = expand(rules.serotyping_seqsero2_wildcards.output.TXT, mode='Allele') if 'fasta' not in config['input'] else [],
        TXT_seqsero2_kmerread = expand(rules.serotyping_seqsero2_wildcards.output.TXT, mode='Kmerread') if 'fasta' not in config['input'] else [],
        JSON_sistr = rules.serotyping_sistr.output.JSON,
        VAL_TSV = rules.serotyping_dump_summary_info.output.VAL_TSV,
        INFORMS_serotyping_sistr = rules.serotyping_sistr.output.INFORMS,
        INFORMS_serotyping_seqsero2 = expand(rules.serotyping_seqsero2_wildcards.output.INFORMS, mode='Kmer')
    output:
        VAL_HTML = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEROTYPE_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'serotyping',
        db_path_sistr = config['serotyping']['sistr']['path'],
        db_path_seqsero2 = config['serotyping']['seqsero2']['path']
    run:
        reportertool = SerovarSalmonellaReporter(camel)
        reportertool.add_input_files({'DIR_sistr': [ToolIODirectory(Path(str(params.db_path_sistr)))]})
        reportertool.add_input_files({'DIR_seqsero2': [ToolIODirectory(Path(str(params.db_path_seqsero2)))]})
        SnakemakeUtils.add_pickle_inputs(reportertool, input, excluded_keys=['VAL_TSV', 'TXT_seqsero2_allele', 'TXT_seqsero2_kmerread'])
        reportertool.add_input_files({'VAL_TSV': [ToolIOFile(Path(input.VAL_TSV))]})
        if input.TXT_seqsero2_allele != []:
            SnakemakeUtils.add_pickle_input(reportertool, 'TXT_seqsero2_allele', Path(input.TXT_seqsero2_allele))
        if input.TXT_seqsero2_kmerread != []:
            SnakemakeUtils.add_pickle_input(reportertool, 'TXT_seqsero2_kmerread', Path(input.TXT_seqsero2_kmerread))
        step = Step(str(rule), reportertool, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reportertool, output)

rule serotyping_report_empty:
    """
    Creates an empty HTML report for the PointFinder analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEROTYPE_REPORT_EMPTY
    params:
        running_dir = Path(config['working_dir']) / 'serotyping'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section(SerovarSalmonellaReporter.TITLE, Path(output.VAL_HTML))
