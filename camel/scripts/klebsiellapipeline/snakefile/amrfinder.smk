from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import assembly_spades

rule amrfinder_run:
    """
    Runs the AMRFinder tool.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA,
        DIR = config['amrfinder']['db']
    output:
        TSV = Path(config['working_dir']) / 'amrfinder' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'amrfinder' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'amrfinder'
    run:
        from camel.app.tools.amrfinder.amrfinder import AMRFinder
        amrfinder = AMRFinder(Camel.get_instance())
        SnakemakeUtils.add_pickle_input(amrfinder, 'FASTA', Path(input.FASTA))
        amrfinder.add_input_files({'DIR': [ToolIODirectory(Path(input.DIR))]})
        step = Step(str(rule), amrfinder, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(amrfinder, output)

rule amrfinder_reporter:
    """
    Creates the output report for AMRFinder.
    """
    input:
        TSV = rules.amrfinder_run.output.TSV,
        INFORMS_amrfinder = rules.amrfinder_run.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / 'amrfinder' / 'html.io'
    params:
        dir_ = Path(config['working_dir']) / 'amrfinder'
    run:
        from camel.app.tools.amrfinder.amrfinderreporter import AMRFinderReporter
        reporter = AMRFinderReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)
