from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import amrfinder

rule amrfinder_run:
    """
    Runs the AMRFinder tool.
    """
    input:
        FASTA = Path(config['working_dir']) / amrfinder.INPUT_AMRFINDER_FASTA,
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

rule amrfinder_report_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'amrfinder' / 'html-empty.io'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('AMRFinder', Path(output.VAL_HTML))

rule amrfinder_dump_summary_info:
    """
    Dumps the summary information for the ResFinder workflow in tabular format.
    """
    input:
        TSV = rules.amrfinder_run.output.TSV
    output:
        TSV = Path(config['working_dir']) / 'amrfinder' / 'summary_amrfinder.tsv'
    run:
        path_tsv = SnakemakeUtils.load_object(Path(input.TSV))[0].path
        data_amr = pd.read_table(path_tsv)

        # Parse perfect & other hits
        data_amr['is_perfect'] = data_amr.apply(lambda x:
            x['% Coverage of reference sequence'] == 100.0 and x['% Identity to reference sequence'] == 100.0, axis=1)
        hits_perfect = list(data_amr[data_amr['is_perfect']]['Gene symbol'])
        hits_other = list(data_amr[~data_amr['is_perfect']]['Gene symbol'])

        # Write to output file
        with open(output.TSV, 'w') as handle:
            handle.write('amrfinder_perfect\t{}'.format(', '.join(hits_perfect) if len(hits_perfect) > 0 else '-'))
            handle.write('\n')
            handle.write('amrfinder_other\t{}'.format(', '.join(hits_other) if len(hits_other) > 0 else '-'))
            handle.write('\n')
