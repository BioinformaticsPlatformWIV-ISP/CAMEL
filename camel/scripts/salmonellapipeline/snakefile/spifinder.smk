import json
from pathlib import Path
from typing import Any, Dict, List, Union, Tuple

import pandas as pd

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.pipelines.salmonella.spifinder import SPIFinder
from camel.app.tools.pipelines.salmonella.spifinderreporter import SPIFinderReporter
from camel.scripts.salmonellapipeline.snakefile import spifinder


camel = Camel.get_instance()


rule spifinder_fastq_run :
    """
    This rule executes spifinder and get the results
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        JSON = Path(config['working_dir']) / spifinder.OUTPUT_JSON_SPIFINDER_FASTQ,
        INFORMS = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_FASTQ_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'spifinder'/ 'spifinder_fastq',
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE',
        db_path = config['spifinder']['path']
    run:
        spifindertool = SPIFinder(camel)
        spifindertool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path)))]})
        if params.read_type == 'PE':
            spifindertool.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        else:
            spifindertool.add_input_files(SnakePipelineUtils.extracts_fq_input(
                Path(input.IO), key_se='FASTQ', read_type=params.read_type))
        step = Step(str(rule), spifindertool, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spifindertool, output)


rule spifinder_fasta_run:
    """
    This rule executes spifinder with fasta of the assembly to obtain  results
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA
    output:
        JSON = Path(config['working_dir']) / spifinder.OUTPUT_JSON_SPIFINDER_FASTA,
        INFORMS = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_FASTA_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'spifinder' / 'spifinder_fasta' ,
        db_path = config['spifinder']['path']
    run:
        spifindertool = SPIFinder(camel)
        spifindertool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path)))]})
        SnakemakeUtils.add_pickle_input(spifindertool, 'FASTA', Path(input.FASTA))
        step = Step(str(rule), spifindertool, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spifindertool, output)


rule create_output_summary_spifinder:
    """
    This rule creates a summmary output for the hits of spifinder in fastq and fasta mode
    DO NOT CHANGE THE ORDER OF THE INPUT FILES!!
    """
    input:
        JSON_FASTQ = rules.spifinder_fastq_run.output.JSON if 'fasta' not in config['input'] else [],
        JSON_FASTA = rules.spifinder_fasta_run.output.JSON,
        INFORMS_spifinder_fastq = rules.spifinder_fastq_run.output.INFORMS if 'fasta' not in config['input'] else [],
        INFORMS_spifinder_fasta = rules.spifinder_fasta_run.output.INFORMS
    output:
        VAL_TSV = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_SUMMARY,
        TSV_documentation = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_DOC,
        JSON = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_SUMMARY_JSON
    params:
        running_dir = Path(config['working_dir']) / 'spifinder'
    run:
        def spifinder_json_parser(json_file_path: Path, tool_informs: Dict[str, Any], mode: str) -> Tuple[List[List[Union[str, int, float], ...]], Dict[Any]]:
            """
            This function is able to parse the output json files of the spifinder tool and returns more favorable outputs
            for the Tsv and the camel Json for Hera.
            :param json_file_path: Path of the json file to be parsed
            :param tool_informs: tool informs corresponding to the run of which json_file_path was the output
            :param mode: fasta or fastq
            :return: a list of hits to be added in the output tsv, a dictionary of hits and metadata to be added to the output json
            """
            with json_file_path.open('r') as file_handle:
                json_file = json.load(file_handle)
            results = []
            spi = json_file['spifinder']["results"]['Salmonella Pathogenicity Islands']['SPI']
            if spi == "No hit found":
                inter_json_dict = { f"spifinder_{mode}" : {'results': results } }
            else:
                hit_dictionary_list = []
                for hits in spi.keys():
                    json_dict = {}
                    header_part1 = ['SPI', 'identity']
                    if mode == 'fasta':
                        header_part2 = ['contig_name', 'positions_in_contig', 'accession', 'insertion_site', 'category_function']
                    else:  # mode == 'fastq':
                        header_part2 = ['accession', 'insertion_site', 'category_function']
                    results.append([spi[hits][hit_property] for hit_property in header_part1] +
                                    [f"{spi[hits]['HSP_length']}/{spi[hits]['template_length']}"] +
                                    [spi[hits][hit_property] for hit_property in header_part2])
                    for hit_property in header_part1 + header_part2:
                        json_dict[hit_property] = spi[hits][hit_property]
                    json_dict['coverage'] = f"{spi[hits]['HSP_length']}/{spi[hits]['template_length']}"
                    hit_dictionary_list.append(json_dict)
                inter_json_dict = { f"spifinder_{mode}" : { 'results' : hit_dictionary_list } }

            inter_json_dict[f"spifinder_{mode}"]['informs_tools'] = { tool_informs.get('_tool', tool_informs['_name']) : {'_name': tool_informs['_name'], '_version': tool_informs['_version'], '_command': tool_informs['_command'], '_tag': tool_informs['_tag']} }
            inter_json_dict[f"spifinder_{mode}"]['informs_dbs'] = {'last_updated': tool_informs['last_update_date'], 'name': tool_informs['key'], 'title': tool_informs['key']}
            return results, inter_json_dict

        with Path(output.VAL_TSV).open('w') as handle:
            meta_json_dict = {}
            if 'fasta' not in config['input']:
                results_fastq_tsv, results_fastq_json = spifinder_json_parser(SnakemakeUtils.load_object(Path(input.JSON_FASTQ))[0].path,
                    SnakemakeUtils.load_object(Path(input[input.INFORMS_spifinder_fastq])), 'fastq')
                meta_json_dict.update(results_fastq_json)
                handle.write(f"spifinder_fastq\t{results_fastq_tsv}\n")
            results_fasta, results_fasta_json = spifinder_json_parser(SnakemakeUtils.load_object(Path(input.JSON_FASTA))[0].path,
                SnakemakeUtils.load_object(Path(input[input.INFORMS_spifinder_fasta])), 'fasta')
            meta_json_dict.update(results_fasta_json)
            handle.write(f"spifinder_fasta\t{results_fasta}\n")
            with Path(output.JSON).open('w') as handle2:
                handle2.write(json.dumps(meta_json_dict))

        # Generate a tsv which documents the meaning of the function categories in the fasta results
        file = pd.read_csv(config['spifinder']['metadata'], delimiter=';')
        file.to_csv(output.TSV_documentation, sep='\t')

rule create_output_report_spifinder:
    """
    This rule creates a simple output report, combining both spifinder tables in one report
    """
    input:
        JSON_FASTQ = rules.spifinder_fastq_run.output.JSON if 'fasta' not in config['input'] else [],
        JSON_FASTA = rules.spifinder_fasta_run.output.JSON,
        VAL_TSV = rules.create_output_summary_spifinder.output.VAL_TSV,
        TSV_documentation = rules.create_output_summary_spifinder.output.TSV_documentation,
        INFORMS_spifinder_fastq = rules.spifinder_fastq_run.output.INFORMS if 'fasta' not in config['input'] else [],
        INFORMS_spifinder_fasta = rules.spifinder_fasta_run.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'spifinder'
    threads: 8
    run:
        spifinder_reporter = SPIFinderReporter(camel)
        SnakemakeUtils.add_pickle_inputs(spifinder_reporter, input, excluded_keys=['VAL_TSV', 'TSV_documentation', 'JSON_FASTQ', 'INFORMS_spifinder_fastq'])
        spifinder_reporter.add_input_files({'TSV_output': [ToolIOFile(Path(input.VAL_TSV))],
                                      'TSV_documentation': [ToolIOFile(Path(input.TSV_documentation))]})
        if input.JSON_FASTQ:
            SnakemakeUtils.add_pickle_input(spifinder_reporter, 'JSON_FASTQ', Path(input.JSON_FASTQ))
        if input.INFORMS_spifinder_fastq:
            SnakemakeUtils.add_pickle_inputs(spifinder_reporter, input, excluded_keys=['VAL_TSV','TSV_documentation', 'JSON_FASTQ', 'JSON_FASTA', 'INFORMS_spifinder_fasta'])
        step = Step(str(rule), spifinder_reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spifinder_reporter, output)


rule spifinder_report_empty:
    """
    Creates an empty HTML report for the PointFinder analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_REPORT_EMPTY
    params:
        running_dir = Path(config['working_dir']) / 'spifinder'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section(SPIFinderReporter.TITLE, Path(output.VAL_HTML))
