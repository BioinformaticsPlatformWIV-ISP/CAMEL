from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.bacilluspipeline.snakefile import ani

camel = Camel.get_instance()

rule fastani_run:
    """
    Runs fastani on assembled contigs.
    """
    input:
        FASTA_Q = Path(config['working_dir']) / ani.INPUT_FASTA_ANI
    output:
        TSV = Path(config['working_dir']) / ani.OUTPUT_VAL_ANI,
        INFORMS = Path(config['working_dir']) / ani.OUTPUT_INFORMS_ANI
    params:
        running_dir = Path(config['working_dir']) / 'ani'
    run:
        from camel.app.tools.fastani.fastani import FastANI
        fastani = FastANI(camel)
        SnakemakeUtils.add_pickle_inputs(fastani, input)
        fastani.add_input_files({'TSV_FASTA_R': [ToolIOFile(Path(config['fastani']['path']))]})
        step = Step(str(rule), fastani, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastani, output)

rule fastani_report:
    """
    Creates the report section for fastani.
    """
    input:
        TSV = rules.fastani_run.output.TSV,
        INFORMS_fastani = rules.fastani_run.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / ani.OUTPUT_ANI_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'ani',
        sample_name = config['sample_name'],
        species = config['species']
    run:
        from camel.app.tools.fastani.fastanireporter import FastANIReporter
        ani_report = FastANIReporter(camel)
        SnakemakeUtils.add_pickle_inputs(ani_report, input)
        ani_report.update_parameters(sample_name=params.sample_name, species=params.species)
        step = Step(str(rule), ani_report, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(ani_report, output)

rule fastani_report_empty:
    """
    Creates an empty HTML report for the FastANI analysis.
    """
    output:
        HTML = Path(config['working_dir']) / ani.OUTPUT_ANI_REPORT_EMPTY
    params:
        running_dir = Path(config['working_dir']) / 'ani'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.fastani.fastanireporter import FastANIReporter
        SnakePipelineUtils.create_empty_report_section(FastANIReporter.TITLE, Path(output.HTML))

rule fastani_dump_summary_info:
    """
    Dumps the summary information for the FastANI workflow in tabular format.
    """
    input:
        TSV = rules.fastani_run.output.TSV
    output:
        TSV = Path(config['working_dir']) / ani.OUTPUT_ANI_SUMMARY
    run:
        tsv_fastani = SnakemakeUtils.load_object(Path(input.TSV))[0].path
        fastani_table = pd.read_table(tsv_fastani, header=None)
        pd.set_option('display.max_columns',None)
        with open(output.TSV, 'w') as handle:
            handle.write('fastani_closest_species\t{}'.format([Path(fastani_table[1][0]).stem, str(fastani_table[2][0])+'%']))
            handle.write('\n')
