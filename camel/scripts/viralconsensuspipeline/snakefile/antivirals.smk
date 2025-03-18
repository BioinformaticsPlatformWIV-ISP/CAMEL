import itertools
import json

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.viralconsensuspipeline.snakefile import antivirals
from camel.scripts.viralconsensuspipeline.snakefile import nextclade3

rule antivirals_check_mutations:
    """
    Queries the detected antiviral mutations against the database.
    """
    input:
        DB = '/db/pipelines/viral_consensus/version_1.1/antivirals',
        TSV = lambda wildcards: nextclade3.get_nextclade_output(wildcards, checkpoints, 'TSV', config),
        INFORMS_subtype = lambda wildcards: nextclade3.get_informs_subtype(wildcards, checkpoints)
    output:
        JSON = Path(config['working_dir']) / 'antivirals' / 'json.io',
    params:
        dir_ = Path(config['working_dir']) / 'antivirals',
        species = config.get('antivirals').get('species')
    run:
        from camel.app.tools.pipelines.viral_consensus.antiviralsdetection import AntiviralsDetection
        detection = AntiviralsDetection(Camel.get_instance())
        detection.add_input_files({
            'DB': [ToolIODirectory(Path(input.DB))],
            'TSV': list(itertools.chain(*[SnakemakeUtils.load_object(Path(tsv)) for tsv in input.TSV])),
        })
        subtype = SnakemakeUtils.load_object(Path(str(input.INFORMS_subtype)))['subtype']
        if subtype is None:
            raise RuntimeError(f'A subtype should be detected to detected antiviral mutations.')
        detection.update_parameters(species=params.species, subtype=subtype)
        detection.run(params.dir_)
        SnakemakeUtils.dump_tool_outputs(detection, output)

rule antivirals_report:
    """
    Creates a HTML report for the antiviral resistance screening.
    """
    input:
        JSON = rules.antivirals_check_mutations.output.JSON
    output:
        VAL_HTML = Path(config['working_dir']) / antivirals.OUTPUT_REPORT
    params:
        dir_ = Path(config['working_dir']) / 'antivirals'
    run:
        from camel.app.tools.pipelines.viral_consensus.antiviralsreporter import AntiviralsReporter
        reporter = AntiviralsReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        reporter.run(params.dir_)
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule antivirals_report_empty:
    """
    Creates an empty output report when the assay is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / antivirals.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Antiviral resistance', Path(output.VAL_HTML))

rule antivirals_summary:
    """
    Creates the summary output for the antiviral resistance screening.
    """
    input:
        JSON = rules.antivirals_check_mutations.output.JSON
    output:
        TSV = Path(config['working_dir']) / 'antivirals' / 'summary.tsv'
    run:
        path_json = SnakemakeUtils.load_object(Path(input.JSON))[0].path
        with open(path_json) as handle:
            data_mutations = json.load(handle)

        data_summary = [
            ('antivirals_mutations', json.dumps(data_mutations['mutations'])),
            ('antivirals_associations', json.dumps(data_mutations['associations'])),
        ]
        with open(output.TSV, 'w') as handle:
            for key, value in data_summary:
                handle.write('\t'.join([key, str(value)]))
                handle.write('\n')
