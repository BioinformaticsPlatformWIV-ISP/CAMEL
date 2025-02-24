from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import assembly
from camel.scripts.mycobacteriumpipeline.snakefile import spoligotyping


rule spoligotyping_downsample:
    """
    Optional downsampling if the coverage is too high.
    """
    input:
        IO_FASTQ = Path(config['working_dir']) / 'fq_dict.io',
        INFORMS_coverage = Path(config['working_dir']) / assembly.get_depth_inform('fastq_pe', 'ref')
    output:
        FASTQ_PE = Path(config['working_dir']) / 'spoligotyping' / 'downsampling' / 'fastq-ds.io',
        INFORMS_spoligo_param = Path(config['working_dir']) / 'spoligotyping' / 'downsampling' / 'informs-param.io'
    params:
        dir_ = Path(config['working_dir']) / 'spoligotyping' / 'downsampling',
        read_key = 'SE' if config['input_type'] == 'ont' else 'PE'
    run:
        import logging
        from camel.app.tools.seqtk.seqtksubsample import SeqtkSubsample
        estimated_coverage = SnakemakeUtils.load_object(Path(input.INFORMS_coverage))['median_depth']
        downsample_factor = 50 / max(estimated_coverage, 1)
        logging.info("Downsampling factor for SpoTyping: {:.2f}".format(downsample_factor))
        if downsample_factor < 1:
            seqtk = SeqtkSubsample(Camel.get_instance())
            seqtk.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO_FASTQ)))
            step = Step(str(rule), seqtk, Camel.get_instance(), Path(params.dir_))
            seqtk.update_parameters(fraction=downsample_factor)
            step.run_step()
            SnakemakeUtils.dump_object({
                'min_strict': 5,
                'min_relaxed': 5,
                'downsample_factor': '{:.2f}'.format(downsample_factor)
            }, Path(output.INFORMS_spoligo_param))
            SnakemakeUtils.dump_object(seqtk.tool_outputs['FASTQ_PE'], Path(output.FASTQ_PE))
        else:
            fq_dict = SnakemakeUtils.load_object(Path(input.IO_FASTQ))
            SnakemakeUtils.dump_object(fq_dict[params.read_key], Path(output.FASTQ_PE))
            SnakemakeUtils.dump_object({
                'min_strict': max(round(estimated_coverage / 10), 3),
                'min_relaxed': max(round(estimated_coverage / 10), 3),
                'downsample_factor': 'NA'
            }, Path(output.INFORMS_spoligo_param))

rule spoligotyping_spotyping:
    """
    Runs the SpoTyping tool.
    """
    input:
        FASTQ = rules.spoligotyping_downsample.output.FASTQ_PE if config['input_type'] in ( 'illumina', 'ont') else [],
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA if config['input_type'] in ('fasta', 'fasta_with_vcf') else [],
        INFORMS_spoligo_param = rules.spoligotyping_downsample.output.INFORMS_spoligo_param if config['input_type'] in ('illumina', 'ont') else []
    output:
        VAL_type_binary = Path(config['working_dir']) / 'spoligotyping' / 'VAL_binary.io',
        VAL_type_octal = Path(config['working_dir']) / 'spoligotyping' / 'VAL_octal.io',
        LOG = Path(config['working_dir']) / 'spoligotyping' / 'log.io',
        INFORMS = Path(config['working_dir']) / spoligotyping.OUTPUT_SPOLIGOTYPING_INFORMS
    params:
        dir_ = Path(config['working_dir']) / 'spoligotyping',
        key = 'FASTQ' if config['input_type'] in ('illumina', 'ont') else 'FASTA'
    run:
        from camel.app.tools.spotyping.spotyping import SpoTyping
        spotyping = SpoTyping(Camel.get_instance())
        if params.key == 'FASTQ':
            SnakemakeUtils.add_pickle_input(spotyping, params.key, Path(input.FASTQ))
            spotyping_params = SnakemakeUtils.load_object(Path(input.INFORMS_spoligo_param))
            spotyping.update_parameters(
                swift='off', min_strict=spotyping_params['min_strict'], min_relaxed=spotyping_params['min_relaxed'])
        else:
            SnakemakeUtils.add_pickle_input(spotyping, params.key, Path(input.FASTA))
            spotyping.update_parameters(swift='off', fasta=None)
        step = Step(str(rule), spotyping, Camel.get_instance(), Path(params.dir_))
        step.run_step()
        spotyping.informs['_tag'] = 'Spoligotyping'
        SnakemakeUtils.dump_tool_outputs(spotyping, output)

rule spoligotyping_report:
    """
    Creates a report for the spoligotyping.
    """
    input:
        VAL_type_binary = rules.spoligotyping_spotyping.output.VAL_type_binary,
        VAL_type_octal = rules.spoligotyping_spotyping.output.VAL_type_octal,
        LOG = rules.spoligotyping_spotyping.output.LOG,
        INFORMS_spotyping = rules.spoligotyping_spotyping.output.INFORMS,
        INFORMS_spoligo_param = rules.spoligotyping_downsample.output.INFORMS_spoligo_param if config['input_type'] not in ('fasta', 'fasta_with_vcf') else []
    output:
        VAL_HTML = Path(config['working_dir']) / spoligotyping.OUTPUT_SPOLIGOTYPING_REPORT,
        INFORMS = Path(config['working_dir']) / 'spoligotyping' / 'informs-report.io'
    params:
        dir_ = Path(config['working_dir'], 'spoligotyping'),
        input_type = config['input_type']
    run:
        from camel.app.tools.spotyping.spotypingreporter import SpoTypingReporter
        reporter = SpoTypingReporter(Camel.get_instance())
        if params.input_type not in ('fasta', 'fasta_with_vcf'):
            SnakemakeUtils.add_pickle_inputs(reporter, input)
        else:
            keys = [k for k in input.keys() if k != 'INFORMS_spoligo_param']
            SnakemakeUtils.add_pickle_inputs(reporter, input, keys=keys)
        step = Step(str(rule), reporter, Camel.get_instance(), Path(params.dir_))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule spoligotyping_report_empty:
    """
    Creates an empty report for the spoligotyping.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / spoligotyping.OUTPUT_SPOLIGOTYPING_REPORT_EMPTY
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.tools.spotyping.spotypingreporter import SpoTypingReporter
        section = SpoTypingReporter.generate_empty_section()
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.VAL_HTML))

rule spoligotyping_dump_summary_info:
    """
    Dumps the summary information from the spoligotyping assay.
    """
    input:
        VAL_type_binary = rules.spoligotyping_spotyping.output.VAL_type_binary,
        VAL_type_octal = rules.spoligotyping_spotyping.output.VAL_type_octal,
        INFORMS = rules.spoligotyping_spotyping.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / spoligotyping.OUTPUT_SPOLIGOTYPING_SUMMARY
    run:
        summary_data = [
            ('spoligotype_binary', SnakemakeUtils.load_object(Path(input.VAL_type_binary))[0].value),
            ('spoligotype_octal', SnakemakeUtils.load_object(Path(input.VAL_type_octal))[0].value),
            ('sit_number', SnakemakeUtils.load_object(Path(input.INFORMS)).get('metadata', {}).get('SIT', 'NA'))
        ]
        with open(output.TSV, 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{str(key)}\t{value}')
                handle.write('\n')
