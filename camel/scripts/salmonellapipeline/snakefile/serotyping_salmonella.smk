from pathlib import Path
import json
import re

from camel.app.camel import Camel
from camel.app.tools.pipelines.salmonella.sistr import Sistr
from camel.app.tools.pipelines.salmonella.seqsero2 import SeqSero2
from camel.app.tools.pipelines.salmonella.seqsero2kmerread import SeqSero2KmerRead
from camel.app.tools.pipelines.salmonella.serovarsalmonellareporter import SerovarSalmonellaReporter
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.resources.snakefile import assembly_spades
from camel.scripts.salmonellapipeline.snakefile import serotyping_salmonella

camel = Camel.get_instance()

rule serotyping_sistr:
    """
    This rule executes sistr to obtain serotpying results using cgMLST
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA
    output:
        TSV = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SISTR_SEROTYPE,
        INFORMS = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEROTYPE_SISTR_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'serotyping' / 'serotyping_sistr',
        db_path_sistr= config['serotyping']['sistr']['path']
    threads: 8
    run:
        sistrtool = Sistr(camel)
        sistrtool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path_sistr)))]})
        SnakemakeUtils.add_pickle_input(sistrtool,'FASTA',Path(input.FASTA))
        step = Step(str(rule),sistrtool,camel,params.running_dir,config)
        sistrtool.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sistrtool,output)

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

rule serotyping_seqsero2_kmer:
    """
    This rule executes SeqSero2 to obtain serotpying results
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA
    output:
        TXT = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEQSERO2_SEROTYPE_KMER,
        INFORMS = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEQSERO2_SEROTYPE_KMER_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'serotyping' / 'serotyping_seqsero2_kmer',
        db_path_seqsero2 = config['serotyping']['seqsero2']['path']
    threads: 8
    run:
        seqserotoolk = SeqSero2(camel)
        seqserotoolk.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path_seqsero2)))]})
        SnakemakeUtils.add_pickle_input(seqserotoolk,'FASTA',Path(input.FASTA))
        step = Step(str(rule),seqserotoolk,camel,params.running_dir,config)
        seqserotoolk.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqserotoolk,output)

rule serotyping_seqsero2_allele:
    """
    This rule executes SeqSero2 to obtain serotpying results
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        TXT = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEQSERO2_SEROTYPE_ALLELE,
        INFORMS = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEQSERO2_SEROTYPE_ALLELE_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'serotyping' / 'serotyping_seqsero2_allele',
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE',
        db_path_seqsero2 = config['serotyping']['seqsero2']['path']
    threads: 8
    run:
        seqserotoola = SeqSero2(camel)
        seqserotoola.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path_seqsero2)))]})
        if params.read_type == 'PE':
            seqserotoola.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        else:
            seqserotoola.add_input_files(SnakePipelineUtils.extracts_fq_input(
                Path(input.IO), key_se='FASTQ', read_type=params.read_type))
        step = Step(str(rule),seqserotoola,camel,params.running_dir,config)
        seqserotoola.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqserotoola,output)


rule serotyping_seqsero2_kmerread:
    """
    This rule executes SeqSero2 to obtain serotpying results
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        TXT = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEQSERO2_SEROTYPE_KMERREAD,
        INFORMS = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEQSERO2_SEROTYPE_KMERREAD_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'serotyping' / 'serotyping_seqsero2_kmerread',
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE',
        db_path_seqsero2 = config['serotyping']['seqsero2']['path']
    threads: 8
    run:
        seqserotoola = SeqSero2KmerRead(camel)
        seqserotoola.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path_seqsero2)))]})
        if params.read_type == 'PE':
            seqserotoola.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        else:
            seqserotoola.add_input_files(SnakePipelineUtils.extracts_fq_input(
                Path(input.IO), key_se='FASTQ', read_type=params.read_type))
        step = Step(str(rule),seqserotoola,camel,params.running_dir,config)
        seqserotoola.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqserotoola,output)

rule create_output_summary_serotyping:
    """
    This rule creates a simple output report for sistr and seqsero2
    """
    input:
        TSV_SISTR = rules.serotyping_sistr.output.TSV,
        INFORMS_SISTR= Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEROTYPE_SISTR_INFORMS,
        TXTSEQSERO2KMER = rules.serotyping_seqsero2_kmer.output.TXT,
        INFORMS_SEQSERO2KMER = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEQSERO2_SEROTYPE_KMER_INFORMS,
        TXTSEQSERO2ALLELE = rules.serotyping_seqsero2_allele.output.TXT if 'fasta' not in config['input'] else [],
        INFORMS_SEQSERO2ALLELE = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEQSERO2_SEROTYPE_ALLELE_INFORMS if 'fasta' not in config['input'] else [],
        TXTSEQSERO2KMERREAD = rules.serotyping_seqsero2_kmerread.output.TXT if 'fasta' not in config['input'] else [],
        INFORMS_SEQSERO2KMERREAD = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEQSERO2_SEROTYPE_KMERREAD_INFORMS if 'fasta' not in config['input'] else []
    output:
        VAL_TSV = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEROTYPE_SUMMARY,
        JSON = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEROTYPE_SUMMARY_JSON
    params:
        running_dir = Path(config['working_dir']) / 'serotyping'
    threads: 8
    run:
        import copy
        meta_json_dict = {}
        #sistr results parsing
        with open(SnakemakeUtils.load_object(Path(input.TSV_SISTR))[0].path) as json_file:
            handle = json.load(json_file)[0]
            with open(output.VAL_TSV,'w') as handle2:
                header_locus = ['Locus', 'serotype_or_group', '% Identity', 'HSP/Locus length', 'Contig', 'Position in contig']
                if handle['qc_status'] == 'PASS':
                    hitstabledict = {'serotype_antigenic_formula':':'.join([str(handle['o_antigen']),
                                                                            str(handle['h1']),
                                                                            str(handle['h2'])
                                                                            ]),
                                     'serotype_serogroup': handle['serogroup'],
                                     'serotype_concensus': handle['serovar'],
                                     'qc_status' : 'PASS'
                                     }
                    hitstabledict_json = copy.deepcopy(hitstabledict)
                    if str(handle['h1_flic_prediction']['is_missing']) == 'False':
                        hit_properties = ['fliC', handle['h1_flic_prediction']['h1'].replace(',', ';'),
                                                                           format(handle['h1_flic_prediction']['top_result']['pident'], '.2f'),
                                                                            '/'.join([str(handle['h1_flic_prediction']['top_result']['length']),
                                                                                      str(handle['h1_flic_prediction']['top_result']['qlen'])]),
                                                                            handle['h1_flic_prediction']['top_result']['stitle'],
                                                                            '...'.join([str(handle['h1_flic_prediction']['top_result']['sstart']),
                                                                                        str(handle['h1_flic_prediction']['top_result']['send'])])
                                                                            ]
                        hitstabledict['hits_serotype_h1_fliC'] = ','.join(hit_properties)
                        json_dict = {}
                        for hit_property_index, hit_property in enumerate(hit_properties):
                            json_dict[header_locus[hit_property_index]] = hit_property
                        hitstabledict_json['hits_serotype_h1_fliC'] = json_dict
                    else:
                        json_dict = {}
                        for item in header_locus:
                            json_dict[item] = "-"
                        hitstabledict_json['hits_serotype_h1_fliC'] = json_dict

                    if str(handle['h2_fljb_prediction']['is_missing']) == 'False':
                        hit_properties = ['fljB', handle['h2_fljb_prediction']['h2'].replace(',', ';'),
                                                                           format(handle['h2_fljb_prediction']['top_result']['pident'], '.2f'),
                                                                            '/'.join([str(handle['h2_fljb_prediction']['top_result']['length']),
                                                                                      str(handle['h2_fljb_prediction']['top_result']['qlen'])]),
                                                                            handle['h2_fljb_prediction']['top_result']['stitle'],
                                                                            '...'.join([str(handle['h2_fljb_prediction']['top_result']['sstart']),
                                                                                        str(handle['h2_fljb_prediction']['top_result']['send'])])
                                                                            ]
                        hitstabledict['hits_serotype_h2_fljB'] = ','.join(hit_properties)
                        json_dict = {}
                        for hit_property_index, hit_property in enumerate(hit_properties):
                            json_dict[header_locus[hit_property_index]] = hit_property
                        hitstabledict_json['hits_serotype_h2_fljB'] = json_dict
                    else:
                        json_dict = {}
                        for item in header_locus:
                            json_dict[item] = "-"
                        hitstabledict_json['hits_serotype_h2_fljB'] = json_dict

                    if str(handle['serogroup_prediction']['wzx_prediction']['is_missing']) == 'False':
                        hit_properties = ['wzx', handle['serogroup_prediction']['wzx_prediction']['serogroup'].replace(',', ';'),
                                                                         format(handle['serogroup_prediction']['wzx_prediction']['top_result']['pident'],'.2f'),
                                                                          '/'.join([str(handle['serogroup_prediction']['wzx_prediction']['top_result']['length']),
                                                                                    str(handle['serogroup_prediction']['wzx_prediction']['top_result']['qlen'])]),
                                                                          handle['serogroup_prediction']['wzx_prediction']['top_result']['stitle'],
                                                                          '...'.join([str(handle['serogroup_prediction']['wzx_prediction']['top_result']['sstart']),
                                                                                      str(handle['serogroup_prediction']['wzx_prediction']['top_result']['send'])]),
                                                                          ]
                        hitstabledict['hits_serotype_o_wzx'] = ','.join(hit_properties)
                        json_dict = {}
                        for hit_property_index, hit_property in enumerate(hit_properties):
                            json_dict[header_locus[hit_property_index]] = hit_property
                        hitstabledict_json['hits_serotype_o_wzx'] = json_dict
                    else:
                        json_dict = {}
                        for item in header_locus:
                            json_dict[item] = "-"
                        hitstabledict_json['hits_serotype_o_wzx'] = json_dict

                    if str(handle['serogroup_prediction']['wzy_prediction']['is_missing']) == 'False':
                        hit_properties = ['wzy', handle['serogroup_prediction']['wzy_prediction']['serogroup'].replace(',', ';'),
                                                                         format(handle['serogroup_prediction']['wzy_prediction']['top_result']['pident'],'.2f'),
                                                                          '/'.join([str(handle['serogroup_prediction']['wzy_prediction']['top_result']['length']),
                                                                                    str(handle['serogroup_prediction']['wzy_prediction']['top_result']['qlen'])]),
                                                                          handle['serogroup_prediction']['wzy_prediction']['top_result']['stitle'],
                                                                          '...'.join([str(handle['serogroup_prediction']['wzy_prediction']['top_result']['sstart']),
                                                                                      str(handle['serogroup_prediction']['wzy_prediction']['top_result']['send'])]),
                                                                          ]
                        hitstabledict['hits_serotype_o_wzy'] = ','.join(hit_properties)
                        json_dict = {}
                        for hit_property_index, hit_property in enumerate(hit_properties):
                            json_dict[header_locus[hit_property_index]] = hit_property
                        hitstabledict_json['hits_serotype_o_wzy'] = json_dict
                    else:
                        json_dict = {}
                        for item in header_locus:
                            json_dict[item] = "-"
                        hitstabledict_json['hits_serotype_o_wzy'] = json_dict

                    handle2.writelines('\t'.join([f"sistr_{k}",v]) + '\n' for k,v in hitstabledict.items())
                else:
                    hitstabledict = {'hits_serotype_h1_fliC': "-",
                                    'hits_serotype_h2_fljB': "-",
                                    'hits_serotype_o_wzx': "-",
                                    'hits_serotype_o_wzy': "-",
                                    'serotype_antigenic_formula': '-',
                                    'serotype_serogroup': '-',
                                    'serotype_concensus': '-',
                                     'qc_status': 'FAIL'
                                     }
                    hitstabledict_json = copy.deepcopy(hitstabledict)
                    for variable in ['hits_serotype_h1_fliC', 'hits_serotype_h2_fljB', 'hits_serotype_o_wzx', 'hits_serotype_o_wzy']:
                        json_dict = {}
                        for item in header_locus:
                            json_dict[item] = "-"
                        hitstabledict_json[variable] = json_dict
                    handle2.writelines('\t'.join([f"sistr_{k}", v]) + '\n' for k, v in hitstabledict.items())
                informs_sistr = SnakemakeUtils.load_object(Path(input.INFORMS_SISTR))
                meta_json_dict.update({'sistr' : {**hitstabledict_json, 'informs_tools' : { informs_sistr['_tool']: {'_name': informs_sistr['_name'], '_version': informs_sistr['_version'], '_command': informs_sistr['_command']}}, 'informs_dbs' : {'last_updated': informs_sistr['last_update_date'], 'name': informs_sistr['key'], 'title': informs_sistr['key']} }})

        #seqsero kmer parsing and allele tool as it's the same file structure
        def seqsero2_parser(seqsero_file:str, seqsero2_mode:str, informs_file: str):
            json_dict = {}
            with open(seqsero_file,'r') as handle:
                seqsero_res = handle.readlines()[2:8]
                seqsero_res = [re.sub(r'([^ ]) ([^ ])',r'\1_\2', res).strip("\n") for res in seqsero_res]
                seqsero_res = [res.replace(':\t','\t') for res in seqsero_res]
                seqsero_res = [seqsero2_mode + "_" + x for x in seqsero_res]
                for res in seqsero_res:
                    json_dict[res.split('\t')[0]] = res.split('\t')[1]
            inter_json_dict = { seqsero2_mode : {**json_dict, 'informs_tools': { informs_file['_tool']: {'_name': informs_file['_name'], '_version': informs_file['_version'], '_command': informs_file['_command'], '_tag': informs_file['_tag']}}, 'informs_dbs' : {'last_updated': informs_file['last_update_date'], 'name': informs_file['key'], 'title': informs_file['key']} }}
            with open(output.VAL_TSV, 'a') as handle:
                handle.writelines(item + '\n' for item in seqsero_res)
            return inter_json_dict

        meta_json_dict.update(seqsero2_parser(SnakemakeUtils.load_object(Path(input.TXTSEQSERO2KMER))[0].path,'seqsero2_kmer', SnakemakeUtils.load_object(Path(input.INFORMS_SEQSERO2KMER))))
        if 'fasta' not in config['input']:
            meta_json_dict.update(seqsero2_parser(SnakemakeUtils.load_object(Path(input.TXTSEQSERO2ALLELE))[0].path,'seqsero2_allele', SnakemakeUtils.load_object(Path(input.INFORMS_SEQSERO2ALLELE))))
            meta_json_dict.update(seqsero2_parser(SnakemakeUtils.load_object(Path(input.TXTSEQSERO2KMERREAD))[0].path,'seqsero2_kmerread', SnakemakeUtils.load_object(Path(input.INFORMS_SEQSERO2KMERREAD))))

        with open(output.JSON, 'w') as handle:
            handle.write(json.dumps(meta_json_dict))

rule create_output_report_serotyping:
    """
    This rule creates a simple output report, combining both serotyping tools
    """
    input:
        TXTSeqSero2kmer = rules.serotyping_seqsero2_kmer.output.TXT,
        TXTSeqSero2allele = rules.serotyping_seqsero2_allele.output.TXT if 'fasta' not in config['input'] else [],
        TXTSeqSero2kmerread = rules.serotyping_seqsero2_kmerread.output.TXT if 'fasta' not in config['input'] else [],
        TSV_SISTR = rules.serotyping_sistr.output.TSV,
        TSV_output = rules.create_output_summary_serotyping.output.VAL_TSV,
        INFORMS_serotyping_sistr = rules.serotyping_sistr.output.INFORMS,
        INFORMS_serotyping_seqsero2 = rules.serotyping_seqsero2_kmer.output.INFORMS
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
        SnakemakeUtils.add_pickle_inputs(reportertool,input,excluded_keys=['TSV_output', 'TXTSeqSero2allele', 'TXTSeqSero2kmerread'])
        reportertool.add_input_files({'TSV_output': [ToolIOFile(Path(input.TSV_output))]})
        if input.TXTSeqSero2allele != []:
            SnakemakeUtils.add_pickle_input(reportertool,'TXTSeqSero2allele',Path(input.TXTSeqSero2allele))
        if input.TXTSeqSero2kmerread != []:
            SnakemakeUtils.add_pickle_input(reportertool,'TXTSeqSero2kmerread',Path(input.TXTSeqSero2kmerread))
        step = Step(str(rule),reportertool,camel,params.running_dir,config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reportertool,output)