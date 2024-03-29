import logging
from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.tools.mash.mashscreen import MashScreen
from camel.scripts.viralconsensuspipeline.snakefile import refselection


rule ref_selection_mash_screen:
    """
    Calculates similarity to genomes in the reference database using mash screen.
    """
    input:
        IO_FASTQ = Path(config['working_dir']) / 'fq_dict.io',
        DB = lambda wildcards: (Path(config['ref_selection']['db'], 'mash', f'{wildcards.segment}.msh')) if config['ref_selection']['db'] is not None else []
    output:
        TSV = Path(config['working_dir']) / 'ref_selection' / 'mash_screen' / '{segment}' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'ref_selection' / 'mash_screen' / '{segment}' / 'informs.io'
    threads: 8
    params:
        dir_working = lambda wildcards: Path(config['working_dir']) / 'ref_selection' / 'mash_screen' / wildcards.segment,
        input_type = config['input_type'],
        segment = lambda wildcards: wildcards.segment
    run:
        from camel.app.io.tooliofile import ToolIOFile

        # Load FASTQ data
        if params.input_type == 'nanopore':
            fq_in = [f.path for f in SnakemakeUtils.load_object(Path(input.IO_FASTQ))['SE']]
        else:
            fq_in = [f.path for f in SnakemakeUtils.load_object(Path(input.IO_FASTQ))['PE']]

        mash = MashScreen(Camel.get_instance())
        mash.add_input_files({
            'FASTQ': [ToolIOFile(p) for p in fq_in],
            'DB': [ToolIOFile(Path(str(input.DB)))]
        })
        mash.update_parameters(threads=threads)
        mash.run(Path(str(params.dir_working)))
        SnakemakeUtils.dump_tool_outputs(mash, output)

rule ref_selection_create_fasta:
    """
    Creates a FASTA file by parsing the mash output and selection the best matching reference for each segment.
    """
    input:
        TSV = ([str(rules.ref_selection_mash_screen.output.TSV).format(segment=seg) for seg in refselection.get_segments(
            Path(config['ref_selection']['db']))]) if config['ref_selection']['db'] is not None else [],
        DB = Path(config['ref_selection']['db']) if config['ref_selection'].get('db') is not None else[]
    output:
        TSV = Path(config['working_dir']) / 'ref_selection' / 'create_fasta' / 'tsv.io',
        JSON = Path(config['working_dir']) / 'ref_selection' / 'create_fasta' / 'json.io',
        FASTA = Path(config['working_dir']) / 'ref_selection' / 'create_fasta' / 'fasta.io'
    params:
        dir_ = Path(config['working_dir']) / 'ref_selection' / 'create_fasta'
    run:
        import itertools
        from camel.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.pipelines.viral_consensus.refselectioncreatefasta import RefSelection
        ref_selection = RefSelection(Camel.get_instance())
        tsv_in = list(itertools.chain(*[SnakemakeUtils.load_object(Path(io)) for io in input.TSV]))
        ref_selection.add_input_files({'TSV': tsv_in, 'DB': [ToolIODirectory(Path(input.DB))]})
        ref_selection.run(params.dir_)
        SnakemakeUtils.dump_tool_outputs(ref_selection, output)

rule ref_selection_report:
    """
    Creates the HTML output report for the reference selection.
    """
    input:
        TSV = rules.ref_selection_create_fasta.output.TSV,
        JSON = rules.ref_selection_create_fasta.output.JSON,
        FASTA = rules.ref_selection_create_fasta.output.FASTA,
        DB = Path(config['ref_selection']['db']) if config['ref_selection'].get('db') is not None else []
    output:
        VAL_HTML = Path(config['working_dir']) / refselection.OUTPUT_REF_SELECTION_REPORT
    params:
        dir_ = Path(config['working_dir']) / 'ref_selection' / 'report'
    run:
        from camel.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.pipelines.viral_consensus.reporterrefselection import ReporterRefSelection
        reporter = ReporterRefSelection(Camel.get_instance())
        step = Step(str(rule), reporter, Camel.get_instance(), Path(params.dir_))
        SnakemakeUtils.add_pickle_inputs(reporter, input, excluded_keys=['DB'])
        reporter.add_input_files({'DB': [ToolIODirectory(Path(input.DB))]})
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule ref_selection_report_empty:
    """
    Creates an empty report when reference genome selection is disabled. 
    """
    output:
        VAL_HTML = Path(config['working_dir']) / refselection.OUTPUT_REF_SELECTION_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Reference selection', Path(output.VAL_HTML))

rule ref_selection_dump_summary_info:
    """
    Creates the summary information for the reference selection.
    """
    input:
        JSON = rules.ref_selection_create_fasta.output.JSON
    output:
        TSV = Path(config['working_dir']) / refselection.OUTPUT_REF_SELECTION_SUMMARY
    run:
        import json

        # Load JSON input
        path_json = SnakemakeUtils.load_object(Path(input.JSON))[0].path
        logging.info(f'Reading: {path_json}')
        with path_json.open() as handle:
            data_ref_selection = json.load(handle)

        # Output the selected references
        rows_out = []
        with open(output.TSV, 'w') as handle:
            for segment, ref in data_ref_selection.items():
                handle.write(f"ref_selection-{segment}\t{ref if ref is not None else '-'}")
                handle.write('\n')
