from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import lrefinder as lrefinder_workflow

rule run_lrefinder:
    """
    Runs the LRE-Finder tool.
    """
    input:
        IO = lrefinder_workflow.get_input(config)
    output:
        INFORMS = Path(config['working_dir']) / 'lrefinder' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'lrefinder',
        input_type = config['input_type']
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.lrefinder.lrefinder import LREFinder
        lrefinder = LREFinder(Camel.get_instance())
        if params.input_type != 'fasta':
            key_reads = 'PE' if params.input_type == 'illumina' else 'SE'
            fq_dict = SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_se = 'FASTQ_SE', drop_empty = True, read_type = key_reads)
            #fq_dict = SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE')
            lrefinder.add_input_files(fq_dict)
        else:
            # When the input type is FASTA the simulated reads are used as input
            SnakemakeUtils.add_pickle_input(lrefinder, 'FASTQ_PE', Path(input.IO))
        step = Step(str(rule), lrefinder, Camel.get_instance(), params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(lrefinder, output)

rule lre_finder_report:
    """
    Creates the HTML output for the LRE-Finder tool.
    """
    input:
        INFORMS_lrefinder = rules.run_lrefinder.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / lrefinder_workflow.OUTPUT_LREFINDER_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'lrefinder',
        input_type = config['input_type']
    run:
        from camel.app.tools.lrefinder.lrefinderreporter import LREFinderReporter
        reporter = LREFinderReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        if params.input_type == 'fasta':
            reporter.update_parameters(pseudo_reads=True)
        step = Step(str(rule), reporter, Camel.get_instance(), params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule lre_finder_report_empty:
    """
    Creates an empty HTML output when the LRE-Finder tool is disabled.
    """
    output:
        HTML = Path(config['working_dir']) / lrefinder_workflow.OUTPUT_LREFINDER_REPORT_EMPTY
    params:
        running_dir = Path(config['working_dir']) / 'lrefinder'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('LRE-Finder', Path(output.HTML))

rule lre_finder_summary:
    """
    Creates the summary output for the LRE-Finder tool.
    """
    input:
        INFORMS = rules.run_lrefinder.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / lrefinder_workflow.OUTPUT_LREFINDER_SUMMARY
    run:
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        with open(output.TSV, 'w') as handle:
            for key in ('species', 'genes', 'mutations'):
                handle.write('\t'.join([f'lrefinder_{key}', str(informs[key])]))
                handle.write('\n')
