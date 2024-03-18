from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import assembly_spades


rule resfinder4_run:
    """
    Runs the ResFinder4 tool.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA,
        DIR = config['resfinder4']['db']
    output:
        TSV_genes = Path(config['working_dir']) / 'resfinder4' / 'tsv-genes.io',
        TSV_point = Path(config['working_dir']) / 'resfinder4' / 'tsv-point.io',
        TSV_pheno_species = Path(config['working_dir']) / 'resfinder4' / 'tsv-pheno-species.io',
        TSV_pheno_general = Path(config['working_dir']) / 'resfinder4' / 'tsv-pheno-general.io',
        INFORMS = Path(config['working_dir']) / 'resfinder4' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'resfinder4',
        species = config['resfinder4'].get('species'),
        point = config['resfinder4'].get('point', True)
    run:
        from camel.app.tools.resfinder.resfinder import ResFinder
        resfinder = ResFinder(Camel.get_instance())
        SnakemakeUtils.add_pickle_input(resfinder, 'FASTA', Path(input.FASTA))
        resfinder.add_input_files({'DIR': [ToolIODirectory(Path(input.DIR))]})
        resfinder.update_parameters(min_cov=0.9, acquired=True, point=params.point)
        if params.species is not None:
            resfinder.update_parameters(species=f'"{params.species}"')
        step = Step(str(rule), resfinder, Camel.get_instance(), params.dir_)
        step.run_step()
        if params.point:
            SnakemakeUtils.dump_tool_outputs(resfinder, output)
        else:
            SnakemakeUtils.dump_tool_outputs(resfinder, output,
                keys=[key for key in output.keys() if key not in ('TSV_point', 'TSV_pheno_species')])
            SnakemakeUtils.dump_object([], Path(output.TSV_point))
            SnakemakeUtils.dump_object([], Path(output.TSV_pheno_species))

rule resfinder4_reporter:
    """
    Creates the output report for ResFinder4.
    """
    input:
        TSV_genes = rules.resfinder4_run.output.TSV_genes,
        TSV_point = rules.resfinder4_run.output.TSV_point,
        TSV_pheno_species = rules.resfinder4_run.output.TSV_pheno_species,
        TSV_pheno_general = rules.resfinder4_run.output.TSV_pheno_general,
        INFORMS_resfinder = rules.resfinder4_run.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / 'resfinder4' / 'html.io'
    params:
        dir_ = Path(config['working_dir']) / 'resfinder4',
        point= config['resfinder4'].get('point',True)
    run:
        from camel.app.tools.resfinder.resfinderreporter import ResFinderReporter
        reporter = ResFinderReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input, excluded_keys = None if params.point else ['TSV_point', 'TSV_pheno_species'])
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule resfinder4_report_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'resfinder4' / 'html-empty.io'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('ResFinder4', Path(output.VAL_HTML))

rule resfinder4_create_summary:
    """
    Creates a tabular summary output for ResFinder4.
    """
    input:
        INFORMS = rules.resfinder4_run.output.INFORMS,
        TSV_genes = rules.resfinder4_run.output.TSV_genes,
        TSV_point = rules.resfinder4_run.output.TSV_point
    output:
        TSV = Path(config['working_dir']) / 'resfinder4' / 'summary_resfinder.tsv'
    params:
        point = config['resfinder4'].get('point',True)
    run:
        tsv_genes = SnakemakeUtils.load_object(Path(input.TSV_genes))[0].path
        if params.point:
            tsv_point = SnakemakeUtils.load_object(Path(input.TSV_point))[0].path
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        with open(output.TSV, 'w') as handle:
            data_genes = pd.read_table(tsv_genes)
            handle.write(f"resfinder4_genes\t{', '.join(list(data_genes['Resistance gene']))}")
            handle.write('\n')
            if params.point:
                data_mutations = pd.read_table(tsv_point)
                handle.write(f"resfinder4_mutations\t{', '.join(list(data_mutations['Mutation']))}")
                handle.write('\n')
            handle.write(f"resfinder4_tool_version\t{informs['_name']}")
            handle.write('\n')
            handle.write(f"resfinder4_db_version\t{informs['db_version_resfinder']}")
            handle.write('\n')
