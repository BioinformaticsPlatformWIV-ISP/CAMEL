import json
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.pipelines.salmonella.seqsero2 import SeqSero2
from camel.app.tools.pipelines.salmonella.seqsero2reporter import SeqSero2Reporter
from camel.app.tools.pipelines.salmonella.sistr import Sistr
from camel.app.tools.pipelines.salmonella.sistrreporter import SistrReporter
from camel.resources.snakefile import assembly
from camel.scripts.salmonellapipeline.snakefile import serotyping_salmonella

camel = Camel.get_instance()

rule serotyping_sistr:
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
        sistrtool = Sistr(camel)
        sistrtool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path_sistr)))]})
        SnakemakeUtils.add_pickle_input(sistrtool, 'FASTA', Path(input.FASTA))
        step = Step(str(rule), sistrtool, camel, Path(str(params.running_dir)), config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sistrtool, output)

rule serotyping_seqsero2_wildcards:
    """
    This rule executes SeqSero2 in kmer mode.
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
        read_key = lambda wildcards: wildcards.input_format, # todo change to include ont
    run:
        seqserotool = SeqSero2(camel)
        seqserotool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path_seqsero2)))]})
        SnakemakeUtils.add_pickle_input(seqserotool, 'FASTA', Path(input.FASTA))
        if params.read_key == 'fastq_pe':
            seqserotool.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        if params.read_key == 'fastq_se':
            seqserotool.add_input_files(SnakePipelineUtils.extracts_fq_input(
                Path(input.IO), key_se='FASTQ', read_type='SE'))
        seqserotool.update_parameters(mode=str(params.mode))
        step = Step(str(rule), seqserotool, camel, Path(str(params.running_dir)), config)
        step.run_step()
        if config['input_type'] == 'hybrid':
            if params.read_key == 'fastq_pe':
                seqserotool.informs['_tag'] = f"{params.mode} - Illumina"
            else:
                seqserotool.informs['_tag'] = f"{params.mode} - ONT"
        SnakemakeUtils.dump_tool_outputs(seqserotool, output)

rule serotyping_dump_summary_info:
    """
    This rule creates a simple output report for SISTR and SeqSero2.
    """
    input:
        JSON_sistr = rules.serotyping_sistr.output.JSON,
        INFORMS_sistr = Path(config['working_dir']) / serotyping_salmonella.OUTPUT_SEROTYPE_SISTR_INFORMS,
        TXT_seqsero2_kmer = lambda wildcards: str(rules.serotyping_seqsero2_wildcards.output.TXT).format(mode='Kmer', input_format=wildcards.input_format),
        INFORMS_seqsero2_kmer = lambda wildcards: str(rules.serotyping_seqsero2_wildcards.output.INFORMS).format(mode='Kmer', input_format=wildcards.input_format),
        TXT_seqsero2_allele = lambda wildcards: str(rules.serotyping_seqsero2_wildcards.output.TXT).format(mode='Allele', input_format=wildcards.input_format) if 'fasta' not in config['input'] else [],
        INFORMS_seqsero2_allele = lambda wildcards: str(rules.serotyping_seqsero2_wildcards.output.INFORMS).format(mode='Allele', input_format=wildcards.input_format) if 'fasta' not in config['input'] else [],
        TXT_seqsero2_kmerread = lambda wildcards: str(rules.serotyping_seqsero2_wildcards.output.TXT).format(mode='Kmerread', input_format=wildcards.input_format) if 'fasta' not in config['input'] else [],
        INFORMS_seqsero2_kmerread = lambda wildcards: str(rules.serotyping_seqsero2_wildcards.output.INFORMS).format(mode='Kmerread', input_format=wildcards.input_format) if 'fasta' not in config['input'] else []
    output:
        VAL_TSV_sistr = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'summary_out_sistr.tsv',
        VAL_TSV_seqsero2 = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'summary_out_seqsero2.tsv'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / wildcards.input_format
    threads: 8
    run:
        import copy

        # parse obligate Sistr output
        with SnakemakeUtils.load_object(Path(str(input.JSON_sistr)))[0].path.open('r') as handle:
            json_data = json.load(handle)[0]
        header_locus = ['Locus', 'serotype_or_group', '% Identity', 'HSP/Locus length', 'Contig', 'Position in contig']
        if json_data['qc_status'] == 'PASS':
            hits_dict_tsv = {'serotype_antigenic_formula':':'.join([str(json_data['o_antigen']),
                                                                    str(json_data['h1']),
                                                                    str(json_data['h2'])
                                                                    ]),
                             'serotype_serogroup': json_data['serogroup'],
                             'serotype_consensus': json_data['serovar'],
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
                            'serotype_consensus': '-',
                             'qc_status': 'FAIL'
                             }
            for variable in ['hits_serotype_h1_fliC', 'hits_serotype_h2_fljB', 'hits_serotype_o_wzx', 'hits_serotype_o_wzy']:
                hits_dict_tsv[variable] = '-'

        informs_sistr = SnakemakeUtils.load_object(Path(str(input.INFORMS_sistr)))
        with Path(output.VAL_TSV_sistr).open('w') as handle:
            for k, v in hits_dict_tsv.items():
                line = f"sistr_{k}\t{v}\n"
                handle.write(line)
            handle.write(f"sistr_tool_version\t{informs_sistr['_name']}\n")
            handle.write(f"sistr_db_version\t{informs_sistr['last_update_date']}\n")

        # parse obligate seqsero2 output
        informs_seqsero2_kmer = SnakemakeUtils.load_object(Path(str(input.INFORMS_seqsero2_kmer)))
        inter_json_dict, tsv_results = serotyping_salmonella.seqsero2_output_parser(SnakemakeUtils.load_object(Path(str(input.TXT_seqsero2_kmer)))[0].path, 'seqsero2_kmer', informs_seqsero2_kmer)
        with Path(output.VAL_TSV_seqsero2).open('w') as handle:
            handle.writelines(item + '\n' for item in tsv_results)

        # parse facultative seqsero2 output
        if 'fasta' not in config['input']:
            for args_tuple in [(SnakemakeUtils.load_object(Path(str(input.TXT_seqsero2_allele)))[0].path, 'seqsero2_allele', SnakemakeUtils.load_object(Path(str(input.INFORMS_seqsero2_allele)))),
                               (SnakemakeUtils.load_object(Path(str(input.TXT_seqsero2_kmerread)))[0].path, 'seqsero2_kmerread', SnakemakeUtils.load_object(Path(str(input.INFORMS_seqsero2_kmerread))))
                                ]:
                with Path(output.VAL_TSV_seqsero2).open('a') as handle:
                    for item in tsv_results:
                        handle.write(item + '\n')

        with Path(output.VAL_TSV_seqsero2).open('a') as handle:
            handle.write(f"seqsero2_tool_version\t{informs_seqsero2_kmer['_name']}\n")
            handle.write(f"seqsero2_db_version\t{informs_seqsero2_kmer['last_update_date']}\n")

rule create_output_report_serotyping_sistr:
    """
    This rule creates a simple output report, combining both serotyping tools
    """
    input:
        JSON_SISTR = rules.serotyping_sistr.output.JSON,
        VAL_TSV = rules.serotyping_dump_summary_info.output.VAL_TSV_sistr,
        INFORMS_serotyping_sistr = rules.serotyping_sistr.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'html_sistr.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / wildcards.input_format,
        db_path_sistr = config['serotyping']['sistr']['path'],
    run:
        reportertool = SistrReporter(camel)
        reportertool.add_input_files({'DIR_sistr': [ToolIODirectory(Path(str(params.db_path_sistr)))]})
        SnakemakeUtils.add_pickle_inputs(reportertool, input, excluded_keys=['VAL_TSV'])
        reportertool.add_input_files({'TSV_output': [ToolIOFile(Path(str(input.VAL_TSV)))]})
        step = Step(str(rule), reportertool, camel, Path(str(params.running_dir)), config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reportertool, output)

rule create_output_report_serotyping_seqsero2:
    """
    This rule creates a simple output report, combining both serotyping tools
    """
    input:
        TXT_seqsero2_kmer = lambda wildcards: str(rules.serotyping_seqsero2_wildcards.output.TXT).format(mode='Kmer', input_format=wildcards.input_format),
        TXT_seqsero2_allele = lambda wildcards: str(rules.serotyping_seqsero2_wildcards.output.TXT).format(mode='Allele', input_format=wildcards.input_format) if 'fasta' not in config['input'] else [],
        TXT_seqsero2_kmerread = lambda wildcards: str(rules.serotyping_seqsero2_wildcards.output.TXT).format(mode='Kmerread', input_format=wildcards.input_format) if 'fasta' not in config['input'] else [],
        VAL_TSV = rules.serotyping_dump_summary_info.output.VAL_TSV_seqsero2,
        INFORMS_serotyping_seqsero2 = lambda wildcards: str(rules.serotyping_seqsero2_wildcards.output.INFORMS).format(mode='Kmer', input_format=wildcards.input_format)
    output:
        VAL_HTML = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'html_seqsero2.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / wildcards.input_format,
        db_path_seqsero2 = config['serotyping']['seqsero2']['path']
    run:
        reportertool = SeqSero2Reporter(camel)
        reportertool.add_input_files({'DIR_seqsero2': [ToolIODirectory(Path(str(params.db_path_seqsero2)))]})
        SnakemakeUtils.add_pickle_inputs(reportertool, input, excluded_keys=['VAL_TSV', 'TXT_seqsero2_allele', 'TXT_seqsero2_kmerread'])
        reportertool.add_input_files({'TSV_output': [ToolIOFile(Path(input.VAL_TSV))]})
        if input.TXT_seqsero2_allele:
            SnakemakeUtils.add_pickle_input(reportertool, 'TXT_seqsero2_allele', Path(str(input.TXT_seqsero2_allele)))
        if input.TXT_seqsero2_kmerread:
            SnakemakeUtils.add_pickle_input(reportertool, 'TXT_seqsero2_kmerread', Path(str(input.TXT_seqsero2_kmerread)))
        step = Step(str(rule), reportertool, camel, Path(str(params.running_dir)), config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reportertool, output)

rule serotyping_report_empty_sistr:
    """
    Creates an empty HTML report for the sistr analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'html_sistr-empty.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / wildcards.input_format
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section(SistrReporter.TITLE, Path(output.VAL_HTML))

rule serotyping_report_empty_seqsero2:
    """
    Creates an empty HTML report for the seqsero2 analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'serotyping' / '{input_format}' / 'html_seqsero2-empty.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping' / wildcards.input_format
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section(SeqSero2Reporter.TITLE, Path(output.VAL_HTML))
