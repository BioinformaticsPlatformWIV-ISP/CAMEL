from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.pipelines.salmonella.genotyphi import Genotyphi
from camel.app.tools.pipelines.salmonella.genotyphireporter import GenotyphiReporter
from camel.scripts.salmonellapipeline.snakefile import genotyphi

camel = Camel.get_instance()


rule genotyphi_run:
    """
    This rule executes genotyphi.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        CSV = Path(config['working_dir']) /  genotyphi.OUTPUT_GENOTYPHI,
        INFORMS = Path(config['working_dir']) / genotyphi.OUTPUT_GENOTYPHI_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'genotyphi' ,
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE',
        db_path = config['genotyphi']['path']
    threads: 8
    run:
        genotyphitool = Genotyphi(camel)
        genotyphitool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path)))]})
        if params.read_type == 'PE':
            genotyphitool.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        else:
            genotyphitool.add_input_files(SnakePipelineUtils.extracts_fq_input(
                Path(input.IO), key_se='FASTQ', read_type=params.read_type))
        step = Step(str(rule), genotyphitool, camel, params.running_dir, config)
        genotyphitool.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(genotyphitool, output)

rule create_output_summary_genotyphi:
    input:
        CSV = rules.genotyphi_run.output.CSV,
        INFORMS_genotyphi = rules.genotyphi_run.output.INFORMS
    output:
        VAL_TSV = Path(config['working_dir']) / genotyphi.OUTPUT_GENOTYPHI_SUMMARY
    params:
        running_dir = Path(config['working_dir']) / 'genotyphi'
    run:
        import re

        import numpy as np
        import pandas as pd

        results = pd.read_csv(SnakemakeUtils.load_object(Path(input.CSV))[0].path)
        informs_genotyphi = SnakemakeUtils.load_object(Path(input.INFORMS_genotyphi))
        # replace all nan by dashes
        results.replace(np.nan, '-' ,inplace=True)
        with open(output.VAL_TSV, 'w') as handle:
            for row_index in range(results.shape[0]):  # resistance information
                for column_index in range(2, 5):
                    key = re.sub("[(\[].*?[)\]]", "",f"genotyphi_{results.iloc[row_index, 1]}_{results.columns[column_index]}").replace(' ','')
                    value = results.iloc[row_index, column_index]
                    # remove text between brackets to get rid of useless details in the columns headers
                    handle.write(f"{key}\t{value}\n".replace(' ',''))
            for column_index in range(10, 19):  # lineage information
                key = f"genotyphi_{results.columns[column_index]}"
                value = results.iloc[0, column_index]
                handle.write(f"{key}\t{value}\n")

rule create_output_report_genotyphi:
    """
    This rule creates a simple output report for the Mykrobe tool using the Genotyphi database.
    """
    input:
        CSV = rules.genotyphi_run.output.CSV,
        VAL_TSV = rules.create_output_summary_genotyphi.output.VAL_TSV,
        INFORMS_genotyphi = rules.genotyphi_run.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / genotyphi.OUTPUT_GENOTYPHI_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'genotyphi'
    run:
        reportertool = GenotyphiReporter(camel)
        SnakemakeUtils.add_pickle_inputs(reportertool, input, excluded_keys=['VAL_TSV'])
        reportertool.add_input_files({'TSV_output': [ToolIOFile(Path(input.VAL_TSV))]})
        step = Step(str(rule), reportertool, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reportertool, output)

rule genotyphi_report_empty:
    """
    Creates an empty HTML report for the PointFinder analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / genotyphi.OUTPUT_GENOTYPHI_REPORT_EMPTY
    params:
        running_dir = Path(config['working_dir']) / 'genotyphi'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section(GenotyphiReporter.TITLE, Path(output.VAL_HTML))
