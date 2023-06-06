from pathlib import Path

from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.bacilluspipeline.snakefile import ani

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
        step = Step(rule,fastani,camel,params.running_dir,config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastani,output)

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
        sample_name= config['sample_name']
    run:
        from camel.app.tools.fastani.fastanireporter import FastANIReporter
        ani_report = FastANIReporter(camel)
        SnakemakeUtils.add_pickle_inputs(ani_report,input)
        ani_report.update_parameters(sample_name=params.sample_name)
        step = Step(rule,ani_report,camel,params.running_dir,config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(ani_report,output)

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
        SnakePipelineUtils.create_empty_report_section(FastANIReporter.TITLE,Path(output.HTML))

rule fastani_dump_summary_info:
    """
    Dumps the summary information for the FastANI workflow in tabular format.
    """
    input:
        INFORMS = Path(config['working_dir']) / rules.fastani_run.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / ani.OUTPUT_ANI_SUMMARY
    run:
        import json
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        data = []
        with open(output.TSV, 'w') as handle:
            handle.write('{}\t{}\n'.format('fastani_typing_results', json.dumps(data)))
