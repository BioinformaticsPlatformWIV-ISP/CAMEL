from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.bacilluspipeline.snakefile import straingst

camel = Camel.get_instance()

rule straingst_kmerize:
    """
    Runs straingst kmerize on reads.
    """
    input:
        FASTQ = Path(config['working_dir']) / 'fq_dict.io'
    output:
        HDF5 = Path(config['working_dir']) / straingst.OUTPUT_HDF5_STRAINGST
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'straingst' / wildcards.read_type,
        read_type = lambda wildcards: wildcards.read_type
    run:
        from camel.app.tools.strainge.straingstkmerize import StrainGSTKmerize
        straingst_kmerize_ = StrainGSTKmerize(camel)
        fq_in = FastqInput.from_fq_dict(Path(input.FASTQ), f'{params.read_type}')
        if params.read_type == 'ont':
            straingst_kmerize_.add_input_files({'FASTQ': [fq_in.se[0]]})
        if params.read_type == 'illumina':
            straingst_kmerize_.add_input_files({'FASTQ': [fq_in.pe[0], fq_in.pe[1]]})
        step = Step(str(rule), straingst_kmerize_, camel, Path(str(params.running_dir)), config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(straingst_kmerize_, output)

rule straingst_run:
    """
    Runs straingst run on the generated HDF5 against the selected database.
    """
    input:
        HDF5 = rules.straingst_kmerize.output.HDF5
    output:
        TSV_STATS = Path(config['working_dir']) / straingst.OUTPUT_STRAINGST_STATS,
        TSV_STRAINS = Path(config['working_dir']) / straingst.OUTPUT_STRAINGST_STRAINS,
        INFORMS = Path(config['working_dir']) / straingst.OUTPUT_INFORMS_STRAINGST
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'straingst' / wildcards.read_type,
        straingst_db = config['straingst']['db']
    run:
        from camel.app.tools.strainge.straingstrun import StrainGSTRun
        straingst_run_ = StrainGSTRun(camel)
        SnakemakeUtils.add_pickle_inputs(straingst_run_, input)
        straingst_run_.add_input_files({'DB_HDF5': [ToolIOFile(Path(params.straingst_db))]})
        step = Step(str(rule), straingst_run_, camel, Path(str(params.running_dir)), config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(straingst_run_, output)

rule straingst_report:
    """
    Creates the report section for strainGST.
    """
    input:
        TSV = rules.straingst_run.output.TSV_STRAINS,
        INFORMS_straingst = rules.straingst_run.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / straingst.OUTPUT_STRAINGST_REPORT
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'straingst' / wildcards.read_type,
        sample_name = config['sample_name'],
        read_type = lambda wildcards: wildcards.read_type
    run:
        from camel.app.tools.strainge.straingstreporter import StrainGSTReporter
        straingst_reporter = StrainGSTReporter(camel)
        SnakemakeUtils.add_pickle_inputs(straingst_reporter, input)
        straingst_reporter.update_parameters(sample_name=params.sample_name, suffix=str(params.read_type))
        step = Step(str(rule), straingst_reporter, camel, Path(str(params.running_dir)), config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(straingst_reporter, output)

rule straingst_report_empty:
    """
    Creates an empty HTML report for the StrainGST analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / straingst.OUTPUT_STRAINGST_REPORT_EMPTY
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'straingst' / wildcards.read_type
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.strainge.straingstreporter import StrainGSTReporter
        SnakePipelineUtils.create_empty_report_section(StrainGSTReporter.TITLE, Path(output.VAL_HTML))

rule straingst_dump_summary_info:
    """
    Dumps the summary information for the StrainGST workflow in tabular format.
    """
    input:
        TSV = rules.straingst_run.output.TSV_STRAINS
    output:
        TSV = Path(config['working_dir']) / straingst.OUTPUT_STRAINGST_SUMMARY
    params:
        read_type = lambda wildcards: wildcards.read_type
    run:
        tsv_straingst = SnakemakeUtils.load_object(Path(input.TSV))[0].path
        straingst_table = pd.read_table(tsv_straingst)
        with open(output.TSV, 'w') as handle:
            handle.write('straingst_closest_strain_{}\t{}\n'.format(params.read_type, straingst_table.iloc[0]["strain"]))
            handle.write('straingst_breadth_of_coverage_{}\t{}\n'.format(params.read_type, straingst_table.iloc[0]["cov"]))
            handle.write('straingst_evenness_{}\t{}\n'.format(params.read_type, straingst_table.iloc[0]["even"]))
            handle.write('straingst_relative_abundance_{}\t{}\n'.format(params.read_type, straingst_table.iloc[0]["rapct"]))
            handle.write('straingst_score_{}\t{}\n'.format(params.read_type, straingst_table.iloc[0]["score"]))
