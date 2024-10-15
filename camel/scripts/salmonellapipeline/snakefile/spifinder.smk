from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import assembly
from camel.scripts.salmonellapipeline.snakefile import spifinder


camel = Camel.get_instance()


rule spifinder_fastq_run :
    """
    This rule executes SPIFinder and returns the results.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        JSON = Path(config['working_dir']) / 'spifinder' / 'spifinder_fastq' / '{input_format}' / 'spifinder_output.io',
        INFORMS = Path(config['working_dir']) / 'spifinder' / 'spifinder_fastq' / '{input_format}' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'spifinder'/ 'spifinder_fastq' / wildcards.input_format,
        read_key = lambda wildcards: wildcards.input_format,
        db_path = config['spifinder']['path']
    run:
        from camel.app.tools.pipelines.salmonella.spifinder import SPIFinder

        spifindertool = SPIFinder(camel)
        spifindertool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path)))]})
        if params.read_key == 'fastq_pe':
            spifindertool.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        else:
            spifindertool.add_input_files(SnakePipelineUtils.extracts_fq_input(
                Path(input.IO), key_se='FASTQ', read_type='SE'))
        step = Step(str(rule), spifindertool, camel, Path(str(params.running_dir)))
        step.run_step()
        if config['input_type'] == 'hybrid':
            if params.read_key == 'fastq_pe':
                spifindertool.informs['_tag'] = 'Illumina'
            else:
                spifindertool.informs['_tag'] = 'ONT'
        SnakemakeUtils.dump_tool_outputs(spifindertool, output)


rule spifinder_fasta_run:
    """
    This rule executes SPIFinder with FASTA input.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA
    output:
        JSON = Path(config['working_dir']) / spifinder.OUTPUT_JSON_SPIFINDER_FASTA,
        INFORMS = Path(config['working_dir']) / spifinder.OUTPUT_SPIFINDER_FASTA_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'spifinder' / 'spifinder_fasta',
        db_path = config['spifinder']['path']
    run:
        from camel.app.tools.pipelines.salmonella.spifinder import SPIFinder

        spifindertool = SPIFinder(camel)
        spifindertool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path)))]})
        SnakemakeUtils.add_pickle_input(spifindertool, 'FASTA', Path(input.FASTA))
        step = Step(str(rule), spifindertool, camel, params.running_dir)
        step.run_step()
        spifindertool.informs['_tag'] = 'FASTA'
        SnakemakeUtils.dump_tool_outputs(spifindertool, output)

rule create_output_summary_spifinder:
    """
    This rule creates a summary output for the hits of SPIFinder in fastq and fasta mode.
    """
    input:
        JSON_FASTQ = rules.spifinder_fastq_run.output.JSON if 'fasta' not in config['input'] else [],
        JSON_FASTA = rules.spifinder_fasta_run.output.JSON,
        INFORMS_spifinder_fastq = rules.spifinder_fastq_run.output.INFORMS if 'fasta' not in config['input'] else [],
        INFORMS_spifinder_fasta = rules.spifinder_fasta_run.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / 'spifinder' / '{input_format}' / 'summary_out.tsv',
        TSV_documentation = Path(config['working_dir']) / 'spifinder' / '{input_format}' / 'spifinder_function_category.tsv'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'spifinder' / wildcards.input_format
    run:
        informs_fasta = SnakemakeUtils.load_object(Path(input.INFORMS_spifinder_fasta))
        with Path(output.TSV).open('w') as handle:
            if 'fasta' not in config['input']:
                informs_fastq = SnakemakeUtils.load_object(Path(input.INFORMS_spifinder_fastq))
                results_fastq_tsv, results_fastq_json = spifinder.spifinder_json_parser(SnakemakeUtils.load_object(Path(input.JSON_FASTQ))[0].path,
                    informs_fastq, 'fastq')
                handle.write(f"spifinder_fastq\t{results_fastq_tsv}\n")
            results_fasta, results_fasta_json = spifinder.spifinder_json_parser(SnakemakeUtils.load_object(Path(input.JSON_FASTA))[0].path,
                informs_fasta, 'fasta')
            handle.write(f"spifinder_fasta\t{results_fasta}\n")
            handle.write(f"spifinder_tool_version\t{informs_fasta['_name']}\n")
            handle.write(f"spifinder_db_version\t{informs_fasta['last_update_date']}\n")

        # Generate a tsv which documents the meaning of the function categories in the fasta results
        file = pd.read_csv(config['spifinder']['metadata'], delimiter=';')
        file.to_csv(output.TSV_documentation, sep='\t')

rule create_output_report_spifinder:
    """
    This rule creates a simple output report, combining both SPIFinder tables in one report.
    """
    input:
        JSON_FASTQ = rules.spifinder_fastq_run.output.JSON if 'fasta' not in config['input'] else [],
        JSON_FASTA = rules.spifinder_fasta_run.output.JSON,
        TSV = rules.create_output_summary_spifinder.output.TSV,
        TSV_documentation = rules.create_output_summary_spifinder.output.TSV_documentation,
        INFORMS_spifinder_fastq = rules.spifinder_fastq_run.output.INFORMS if 'fasta' not in config['input'] else [],
        INFORMS_spifinder_fasta = rules.spifinder_fasta_run.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / 'spifinder' / '{input_format}' / 'html.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'spifinder' / wildcards.input_format
    threads: 8
    run:
        from camel.app.tools.pipelines.salmonella.spifinderreporter import SPIFinderReporter

        spifinder_reporter = SPIFinderReporter(camel)
        SnakemakeUtils.add_pickle_inputs(spifinder_reporter, input, excluded_keys=['TSV', 'TSV_documentation', 'JSON_FASTQ', 'INFORMS_spifinder_fastq'])
        spifinder_reporter.add_input_files({'TSV_output': [ToolIOFile(Path(input.TSV))],
                                            'TSV_documentation': [ToolIOFile(Path(input.TSV_documentation))]})
        if input.JSON_FASTQ:
            SnakemakeUtils.add_pickle_input(spifinder_reporter, 'JSON_FASTQ', Path(input.JSON_FASTQ))
        if input.INFORMS_spifinder_fastq:
            spifinder_reporter.add_input_informs({'spifinder_fastq': SnakemakeUtils.load_object(Path(str(input.INFORMS_spifinder_fastq)))})
        step = Step(str(rule), spifinder_reporter, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spifinder_reporter, output)

rule spifinder_report_empty:
    """
    Creates an empty HTML report for the SPIFinder analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'spifinder' / '{input_format}' / 'html-empty.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'spifinder' / wildcards.input_format
    run:
        from camel.app.tools.pipelines.salmonella.spifinderreporter import SPIFinderReporter

        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section(SPIFinderReporter.TITLE, Path(output.VAL_HTML))
