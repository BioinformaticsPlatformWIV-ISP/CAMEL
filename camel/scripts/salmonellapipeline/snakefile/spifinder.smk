from pathlib import Path
import pandas as pd
import json
from camel.app.camel import Camel
from camel.app.tools.pipelines.salmonella.spifinder import SPIFinder
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.salmonellapipeline.snakefile import spifinder
from camel.app.tools.pipelines.salmonella.spifinderreporter import SPIFinderReporter


camel = Camel.get_instance()


rule spifinder_fastq_run :
    """
    This rule executes spifinder and get the results
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        JSON = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_FASTQ ,
        INFORMS = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_FASTQ_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'spifinder'/ 'spifinder_fastq' ,
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE' ,
        db_path = config['spifinder']['path']
    run:
        spifinderTool = SPIFinder(camel)
        spifinderTool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path)))]})
        if params.read_type == 'PE':
            spifinderTool.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        else:
            spifinderTool.add_input_files(SnakePipelineUtils.extracts_fq_input(
                Path(input.IO), key_se='FASTQ', read_type=params.read_type))
        step = Step(str(rule),spifinderTool,camel,params.running_dir,config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spifinderTool,output)


rule spifinder_fasta_run:
    """
    This rule executes spifinder with fasta of the assembly to obtain  results
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA
    output:
        JSON = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_FASTA,
        INFORMS = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_FASTA_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'spifinder' / 'spifinder_fasta' ,
        db_path = config['spifinder']['path']
    run:
        spifinderTool = SPIFinder(camel)
        spifinderTool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path)))]})
        SnakemakeUtils.add_pickle_input(spifinderTool,'FASTA',Path(input.FASTA))
        step = Step(str(rule),spifinderTool,camel,params.running_dir,config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spifinderTool,output)


rule create_output_summary_spifinder:
    """
    this rule creates a summmary output for the hits of spifinder in fastq and fasta mode
    DO NOT CHANGE THE ORDER OF THE INPUT FILES!!
    """

    input:
        JSONFASTQ = rules.spifinder_fastq_run.output.JSON if 'fasta' not in config['input'] else [],
        JSONFASTA = rules.spifinder_fasta_run.output.JSON,
        INFORMS_spifinder_fastq = rules.spifinder_fastq_run.output.INFORMS if 'fasta' not in config['input'] else[],
        INFORMS_spifinder_fasta = rules.spifinder_fasta_run.output.INFORMS
    output:
        VAL_TSV = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_SUMMARY,
        DOC_TSV = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_DOC,
        JSON = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_SUMMARY_JSON
    params:
        running_dir = Path(config['working_dir']) / 'spifinder'

    """
    start with the output of the fasta file = json
    """
    run:
        def SPIFinder_parser(filePath, mode):
            with open(filePath) as json_file:
                handle = json.load(json_file)
                results = []
                SPI = handle['spifinder']["results"]['Salmonella Pathogenicity Islands']['SPI']
                if SPI == "No hit found":
                    results = []
                    inter_json_dict = { f"spifinder_{mode}" : {'results': [] } }
                    if mode == 'fasta':
                        input_position = 3
                        if 'fasta' in config['input']:
                            input_position = 1
                    elif mode == 'fastq':
                        input_position = 2
                else:
                    json_dicts_list = []
                    for hits in SPI.keys():
                        json_dict = {}
                        header1 = ['SPI', 'identity']
                        if mode == 'fasta':
                            header2 = ['contig_name', 'positions_in_contig', 'accession', 'insertion_site', 'category_function']
                            input_position = 3
                            if 'fasta' in config['input']:
                                input_position = 1
                        elif mode == 'fastq':
                            header2 = ['accession', 'insertion_site', 'category_function']
                            input_position = 2
                        results.append([SPI[hits][hit_property] for hit_property in header1] +
                                        [f"{SPI[hits]['HSP_length']}/{SPI[hits]['template_length']}"] +
                                        [SPI[hits][hit_property] for hit_property in header2])
                        for hit_property in header1 + header2:
                            json_dict[hit_property] = SPI[hits][hit_property]
                        json_dict['coverage'] = f"{SPI[hits]['HSP_length']}/{SPI[hits]['template_length']}"
                        json_dicts_list.append(json_dict)
                    inter_json_dict = { f"spifinder_{mode}" : { 'results' : json_dicts_list } }

                informs_spifinder = SnakemakeUtils.load_object(Path(input[input_position]))
                inter_json_dict[f"spifinder_{mode}"]['informs_tools'] = { informs_spifinder['_tool'] : {'_name': informs_spifinder['_name'], '_version': informs_spifinder['_version'], '_command': informs_spifinder['_command'], '_tag': informs_spifinder['_tag']} }
                inter_json_dict[f"spifinder_{mode}"]['informs_dbs'] = {'last_updated': informs_spifinder['last_update_date'], 'name': informs_spifinder['key'],'title': informs_spifinder['key']}
            return results, inter_json_dict

                # parse and write results into files
        with open(output.VAL_TSV, 'w') as handle:
            meta_json_dict = {}
            if 'fasta' not in config['input']:
                resultsFASTQ, resultsFASTQ_JSON = SPIFinder_parser(SnakemakeUtils.load_object(Path(input.JSONFASTQ))[0].path,'fastq')
                meta_json_dict.update(resultsFASTQ_JSON)
                handle.write(f"spifinder_fastq\t{resultsFASTQ}\n")
            resultsFASTA, resultsFASTA_JSON = SPIFinder_parser(SnakemakeUtils.load_object(Path(input.JSONFASTA))[0].path,'fasta')
            meta_json_dict.update(resultsFASTA_JSON)
            handle.write(f"spifinder_fasta\t{resultsFASTA}\n")
            with open(output.JSON, 'w') as handle2:
                handle2.write(json.dumps(meta_json_dict))

        # generate a tsv of the function category of the tool to provide in the report
        def TSVmakerDoc(savePath):
            file=pd.read_csv(config['spifinder']['metadata'], delimiter=';')
            file.to_csv(savePath,sep='\t')
        TSVmakerDoc(output.DOC_TSV)


rule create_output_report_spifinder:
    """
    This rule creates a simple output report, combining both spifinder tables in one report
    """
    input:
        JSONFASTQ = rules.spifinder_fastq_run.output.JSON if 'fasta' not in config['input'] else [],
        JSONFASTA = rules.spifinder_fasta_run.output.JSON,
        TSV_output = rules.create_output_summary_spifinder.output.VAL_TSV,
        TSV_doc = rules.create_output_summary_spifinder.output.DOC_TSV,
        INFORMS_spifinder_fastq = rules.spifinder_fastq_run.output.INFORMS if 'fasta' not in config['input'] else [],
        INFORMS_spifinder_fasta = rules.spifinder_fasta_run.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'spifinder',
    threads: 8
    run:
        reportertool = SPIFinderReporter(camel)
        SnakemakeUtils.add_pickle_inputs(reportertool,input,excluded_keys=['TSV_output','TSV_doc', 'JSONFASTQ', 'INFORMS_spifinder_fastq'])
        reportertool.add_input_files({'TSV_output': [ToolIOFile(Path(input.TSV_output))],
                                      'TSV_doc': [ToolIOFile(Path(input.TSV_doc))]})
        if input.JSONFASTQ != []:
            SnakemakeUtils.add_pickle_input(reportertool,'JSONFASTQ',Path(input.JSONFASTQ))
        if input.INFORMS_spifinder_fastq != []:
            SnakemakeUtils.add_pickle_inputs(reportertool,input,excluded_keys=['TSV_output','TSV_doc', 'JSONFASTQ', 'JSONFASTA', 'INFORMS_spifinder_fasta'])
        step = Step(str(rule),reportertool,camel,params.running_dir,config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reportertool,output)


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