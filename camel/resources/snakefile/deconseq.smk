from pathlib import Path

from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import deconseq
from camel.app.components.files.fastqutils import FastqUtils
from camel.app.io.tooliofile import ToolIOFile


rule deconseq_run:
    """
    Reads decontamination using Deconseq.
    """
    input:
        IO = Path(config['working_dir']) / deconseq.INPUT_DECONSEQ_FASTQ
    output:
        FASTQ_PE_CLEAN = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_CLEAN_PE,
        FASTQ_SE_FWD_CLEAN = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_CLEAN_SE_FWD,
        FASTQ_SE_REV_CLEAN= Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_CLEAN_SE_REV,
        INFORMS = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_INFORMS,
        INFORMS_PE_FWD = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_INFORMS_PE_FWD,
        INFORMS_PE_REV = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_INFORMS_PE_REV,
        INFORMS_SE_FWD = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_INFORMS_SE_FWD,
        INFORMS_SE_REV = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_INFORMS_SE_REV,
    params:
        running_dir = Path(config['working_dir']) / 'deconseq',
    threads: 8
    priority: 1
    run:
        from camel.app.tools.deconseq.deconseq import Deconseq
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils


        # Reformat FASTQ dictionary
        fq_dict = SnakePipelineUtils.extracts_fq_input(input.IO, key_pe='FASTQ_PE', key_se='FASTQ_SE')
        fwd_pe_deconseq = Deconseq(camel)
        fwd_pe_deconseq.add_input_files({'FASTQ': [fq_dict['FASTQ_PE'][0]]})
        fwd_pe_deconseq.update_parameters(dbs=config['decontamination']['dbs'])
        if config['decontamination']['dbs_retain']:
            fwd_pe_deconseq.update_parameters(dbs_retain=config['decontamination']['dbs_retain'])
        fwd_pe_step = Step(rule, fwd_pe_deconseq, camel, params.running_dir, config)
        fwd_pe_step.run_step()
        SnakemakeUtils.dump_object(fwd_pe_step.informs, output.INFORMS_PE_FWD)

        rev_pe_deconseq = Deconseq(camel)
        rev_pe_deconseq.add_input_files({'FASTQ': [fq_dict['FASTQ_PE'][1]]})
        rev_pe_deconseq.update_parameters(dbs=config['decontamination']['dbs'])
        if config['decontamination']['dbs_retain']:
            rev_pe_deconseq.update_parameters(dbs_retain=config['decontamination']['dbs_retain'])
        rev_pe_step = Step(rule, rev_pe_deconseq, camel, params.running_dir, config)
        rev_pe_step.run_step()
        SnakemakeUtils.dump_object(rev_pe_step.informs, output.INFORMS_PE_REV)

        if len(fq_dict['FASTQ_SE']) > 0 and '_1U.fastq' in fq_dict['FASTQ_SE'][0].path:
            fwd_se_deconseq = Deconseq(camel)
            fwd_se_deconseq.add_input_files({'FASTQ': [fq_dict['FASTQ_SE'][0]]})
            fwd_se_deconseq.update_parameters(dbs=config['decontamination']['dbs'])
            if config['decontamination']['dbs_retain']:
                fwd_se_deconseq.update_parameters(dbs_retain=config['decontamination']['dbs_retain'])
            fwd_se_step = Step(rule, fwd_se_deconseq, camel, params.running_dir, config)
            fwd_se_step.run_step()
            SnakemakeUtils.dump_object(fwd_se_step.informs, output.INFORMS_SE_FWD)
        else:
            fwd_se_deconseq = None
            SnakemakeUtils.dump_object(None, output.INFORMS_SE_FWD)

        if (len(fq_dict['FASTQ_SE']) > 1 and '_2U.fastq' in fq_dict['FASTQ_SE'][1].path) or (len(fq_dict['FASTQ_SE']) > 0 and '_2U.fastq' in fq_dict['FASTQ_SE'][0].path):
            rev_se_deconseq = Deconseq(camel)
            rev_se_deconseq.add_input_files({'FASTQ': [fq_dict['FASTQ_SE'][1]]})
            rev_se_deconseq.update_parameters(dbs=config['decontamination']['dbs'])
            if config['decontamination']['dbs_retain']:
                rev_se_deconseq.update_parameters(dbs_retain=config['decontamination']['dbs_retain'])
            rev_se_step = Step(rule, rev_se_deconseq, camel, params.running_dir, config)
            rev_se_step.run_step()
            SnakemakeUtils.dump_object(rev_se_step.informs, output.INFORMS_SE_REV)
        else:
            rev_se_deconseq = None
            SnakemakeUtils.dump_object(None, output.INFORMS_SE_REV)

        if config['decontamination']['dbs_retain']:
            fwd_pe = [fwd_pe_deconseq.tool_outputs['FASTQ_Clean'][0].path, fwd_pe_deconseq.tool_outputs['FASTQ_Both'][0].path]
            rev_pe = [rev_pe_deconseq.tool_outputs['FASTQ_Clean'][0].path, rev_pe_deconseq.tool_outputs['FASTQ_Both'][0].path]
            if fwd_se_deconseq:
                se_fwd = [fwd_se_deconseq.tool_outputs['FASTQ_Clean'][0].path, fwd_se_deconseq.tool_outputs['FASTQ_Both'][0].path]
            else:
                se_fwd = []
            if rev_se_deconseq:
                se_rev = [rev_se_deconseq.tool_outputs['FASTQ_Clean'][0].path, rev_se_deconseq.tool_outputs['FASTQ_Both'][0].path]
            else:
                se_rev = []
        else:
            fwd_pe = [fwd_pe_deconseq.tool_outputs['FASTQ_Clean'][0].path]
            rev_pe = [rev_pe_deconseq.tool_outputs['FASTQ_Clean'][0].path]
            fwd_se = [fwd_se_deconseq.tool_outputs['FASTQ_Clean'][0].path] if fwd_se_deconseq else []
            rev_se = [rev_se_deconseq.tool_outputs['FASTQ_Clean'][0].path] if rev_se_deconseq else []

        deconseq_wd = Path(config['working_dir']) / 'deconseq'

        FastqUtils.process_paired_end_se(fwd_pe, rev_pe, fwd_se, rev_se, deconseq_wd / 'deconseq_cleaned_pe_fwd.fq',
                                         deconseq_wd / 'deconseq_cleaned_pe_rev.fq',
                                         deconseq_wd / 'deconseq_cleaned_se_fwd.fq',
                                         deconseq_wd / 'deconseq_cleaned_se_rev.fq')

        SnakemakeUtils.dump_object([ToolIOFile(deconseq_wd / 'deconseq_cleaned_pe_fwd.fq'),
                                    ToolIOFile(deconseq_wd / 'deconseq_cleaned_pe_rev.fq')], output.FASTQ_PE_CLEAN)
        SnakemakeUtils.dump_object([ToolIOFile(deconseq_wd / 'deconseq_cleaned_se_fwd.fq')], output.FASTQ_SE_FWD_CLEAN)
        SnakemakeUtils.dump_object([ToolIOFile(deconseq_wd / 'deconseq_cleaned_se_rev.fq')], output.FASTQ_SE_REV_CLEAN)
        deconseq_informs = deconseq.combine_deconseq_informs(fwd_pe_deconseq, rev_pe_deconseq, fwd_se_deconseq, rev_se_deconseq)
        deconseq_informs['combined'] = {'remaining_pe_reads': FastqUtils.count_reads(str(deconseq_wd / 'deconseq_cleaned_pe_fwd.fq')),
                                        'remaining_se_reads': FastqUtils.count_reads(str(deconseq_wd / 'deconseq_cleaned_se.fq')),
                                        'processed_dbs': deconseq.get_processed_dbs(deconseq_informs)}
        SnakemakeUtils.dump_object(deconseq_informs, deconseq.OUTPUT_DECONSEQ_INFORMS)

rule deconseq_report:
    """
    Creates the HTML report for the decontamination.
    """
    input:
        INFORMS_deconseq = rules.deconseq_run.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'deconseq' / 'report',
        sample_name = config['sample_name']
    run:
        from camel.app.tools.pipelines.deconseq.reporterdeconseq import ReporterDeconseq

        reporter = ReporterDeconseq(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)


rule deconseq_export_summary:
    """
    Exports the summary information for deconseq.
    """
    input:
        INFORMS = rules.deconseq_run.output.INFORMS
    output:
        TSV = Path(config['working_dir'], deconseq.OUTPUT_DECONSEQ_SUMMARY)
    run:

        def get_deconseq_summary_read_count(informs_, deconseq_db_):
            if deconseq_db_ in informs_['deconseq_stats']:
                return informs_['deconseq_stats'][deconseq_db_]['input_reads_count'], informs_['deconseq_stats'][deconseq_db_]['removed_reads_count']
            else:
                return 0, 0

        deconseq_informs = SnakemakeUtils.load_object(input.INFORMS)
        with open(output.TSV, 'w') as handle:
            se_counts = 0
            for read_type in ['PE_FWD', 'PE_REV', 'SE_FWD', 'SE_REV']:
                if deconseq_informs[read_type] is not None:
                    handle.write(f'deconseq_{read_type.lower()}_initial\t{deconseq_informs[read_type]["initial_reads_count"]}\n')
                    handle.write(f'deconseq_{read_type.lower()}_final\t{deconseq_informs[read_type]["final_reads_count"]}\n')
                else:
                    handle.write(f'deconseq_{read_type.lower()}_initial\t0\n')
                    handle.write(f'deconseq_{read_type.lower()}_final\t0\n')
            for deconseq_db in deconseq_informs['combined']['processed_dbs']:
                for read_type in ['PE_FWD', 'PE_REV', 'SE_FWD', 'SE_REV']:
                    if deconseq_informs[read_type] is not None:
                        stats = get_deconseq_summary_read_count(deconseq_informs[read_type], deconseq_db)
                        handle.write(f'deconseq_{read_type.lower()}_{deconseq_db}_input\t{stats[0]}\n')
                        handle.write(f'deconseq_{read_type.lower()}_{deconseq_db}_removed\t{stats[1]}\n')
                    else:
                        handle.write(f'deconseq_{read_type.lower()}_{deconseq_db}_input\t0\n')
                        handle.write(f'deconseq_{read_type.lower()}_{deconseq_db}_removed\t0\n')

rule trimming_deconseq_to_dict:
    """
    Combines the deconseq cleand reads into a dictionary.
    """
    input:
        FASTQ_PE = rules.deconseq_run.output.FASTQ_PE_CLEAN,
        FASTQ_SE_FWD = rules.deconseq_run.output.FASTQ_SE_FWD_CLEAN,
        FASTQ_SE_REV = rules.deconseq_run.output.FASTQ_SE_REV_CLEAN
    output:
        IO = Path(config['working_dir']) / deconseq.OUTPUT_DECONSEQ_DICT
    run:
        output_dict = {
            'PE': SnakemakeUtils.load_object(input.FASTQ_PE)
        }
        se_fwd = SnakemakeUtils.load_object(input.FASTQ_SE_FWD)
        if se_fwd[0].size > 0:
            output_dict['SE_FWD'] = se_fwd
        se_rev = SnakemakeUtils.load_object(input.FASTQ_SE_REV)
        if se_rev[0].size > 0:
            output_dict['SE_REV'] = se_rev
        SnakemakeUtils.dump_object(output_dict, output.IO)