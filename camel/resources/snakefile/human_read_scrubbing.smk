from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import human_read_scrubbing
from camel.resources.snakefile.human_read_scrubbing import get_removed

camel = Camel.get_instance()


rule scrubbing_fasta_fa2fq:
    """
    Convert the input FASTA to FASTQ to be able to be used by the scrubber, it only accepts FASTQ.
    """
    input:
        FASTA = Path(config['working_dir']) / human_read_scrubbing.INPUT_SCRUBBING_FASTA
    output:
        FASTQ_from_fasta = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'input' / 'fasta2fastq.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'human_read_scrubbing' / wildcards.input_format / 'input'
    run:
        from camel.app.components.files.fastautils import FastaUtils

        fasta_path_in = (SnakemakeUtils.load_object(Path(input.FASTA)))[0].path
        fastq_path_out = Path(str(params.running_dir), f"{fasta_path_in.stem}.fastq")

        FastaUtils.convert_fasta_to_fastq(fasta_path_in, fastq_path_out)
        SnakemakeUtils.dump_object([ToolIOFile(fastq_path_out)], Path(output.FASTQ_from_fasta))

rule scrubbing_fastq_interleave_and_gunzip:
    """
    Gunzips input FASTQ files if they are gzipped, and interleaves the input if there are two files, 
    because the tool only accepts a single input FASTQ file.
    """
    input:
        FASTQ = Path(config['working_dir']) / human_read_scrubbing.INPUT_SCRUBBING_FASTQ
    output:
        FASTQ_SINGLE_GUNZIP = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'input' / 'fastq_gunzip_interleaved.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'human_read_scrubbing' / wildcards.input_format / 'input'
    threads: 8
    run:
        from camel.app.components.filesystemhelper import FileSystemHelper
        from camel.app.components.files.fastqutils import FastqUtils

        # Get the FASTQ file(s)
        fastq_in = SnakemakeUtils.load_object(Path(input.FASTQ))
        nb_of_fq_files = len(fastq_in)

        if nb_of_fq_files == 1:
            path_in = fastq_in[0].path
            path_out = Path(str(params.running_dir), path_in.name.replace('.gz', ''))

            if not FileSystemHelper.is_gzipped(path_in):
                SnakemakeUtils.dump_object([ToolIOFile(path_in)], Path(output.FASTQ_SINGLE_GUNZIP))
            else:
                FileSystemHelper.pigz_extract(path_in, path_out, threads=threads)
                SnakemakeUtils.dump_object([ToolIOFile(path_out)], Path(output.FASTQ_SINGLE_GUNZIP))

        else:
            interleaved_out = Path(str(params.running_dir), f"{fastq_in[0].path.stem}_interleaved.fastq")
            FastqUtils.convert_fastqs_to_interleaved_fastq(fastq_in[0].path, fastq_in[1].path, Path(str(params.running_dir), f"{fastq_in[0].path.stem}_interleaved.fastq"))
            SnakemakeUtils.dump_object([ToolIOFile(interleaved_out)], Path(output.FASTQ_SINGLE_GUNZIP))

rule scrubbing_run_scrubber:
    """
    Runs the NCBI human read scrubber on the input FASTQ, only accepts single gunzipped FASTQ files.
    """
    input:
        FASTQ = Path(config['working_dir']) / human_read_scrubbing.INPUT_SCRUBBING_FASTQ if 'fasta' not in config['input'] else [],
        FASTQ_SINGLE_GUNZIP = rules.scrubbing_fastq_interleave_and_gunzip.output.FASTQ_SINGLE_GUNZIP if 'fasta' not in config['input'] else rules.scrubbing_fasta_fa2fq.output.FASTQ_from_fasta
    output:
        FASTQ_SCRUBBED = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'scrubbing' / 'fastq_scrubbed.io',
        FASTQ_REMOVED = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'scrubbing' / 'fastq_removed.io' if config['read_scrubbing'].get('export_removed_reads') else [],
        INFORMS = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'scrubbing' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'human_read_scrubbing' / wildcards.input_format / 'scrubbing',
        input_format = lambda wildcards: wildcards.input_format,
        export_removed_reads = config['read_scrubbing'].get('export_removed_reads')
    threads: max(16, workflow.cores * 0.75)
    run:
        from camel.app.tools.ncbihumanreadscrubber.ncbihumanreadscrubber import NcbiHumanReadScrubber
        if params.input_format != 'fasta':
            # Get the FASTQ file(s)
            fastq_in = SnakemakeUtils.load_object(Path(input.FASTQ))
            fqfile_number = len(fastq_in)
            interleaved = True if fqfile_number == 2 else False
        else:
            interleaved = False

        outputfile_scrubbing = str(Path(output.FASTQ_SCRUBBED).with_suffix('.fastq')) if params.input_format != 'fasta' \
            else str(Path(str(params.running_dir), f"{SnakemakeUtils.load_object(Path(input.FASTQ_SINGLE_GUNZIP))[0].path.name}_scrubbed.fastq"))

        scrubber = NcbiHumanReadScrubber(camel)
        step = Step(str(rule), scrubber, camel, Path(str(params.running_dir)))
        scrubber.update_parameters(interleaved=interleaved, outputfile=outputfile_scrubbing, threads=threads)
        if params.export_removed_reads:
            outputfile_removed_reads = str(Path(output.FASTQ_REMOVED).with_suffix('.fastq'))
            scrubber.update_parameters(export_human_reads=params.export_removed_reads, outputfile_removed=outputfile_removed_reads)
        SnakemakeUtils.add_pickle_inputs(scrubber, input, excluded_keys=['FASTQ'])
        step.run_step()

        if config['input_type'] == 'hybrid':
            if params.input_format == 'fastq_pe':
                scrubber.informs['_tag'] = 'Illumina'
            else:
                scrubber.informs['_tag'] = 'ONT'

        SnakemakeUtils.dump_tool_outputs(scrubber, output, keys=['FASTQ_SCRUBBED', 'INFORMS'])
        SnakemakeUtils.dump_tool_outputs(scrubber, output, keys=['FASTQ_REMOVED'], ignore_missing_output=True) if params.export_removed_reads else ''

rule scrubbing_fasta_fq2fa:
    """
    Convert the FASTQ back to FASTA in order for the rest of the pipeline to be able 
    to use it as output or as input in bacterial pipelines.
    """
    input:
        FASTQ_SCRUBBED = rules.scrubbing_run_scrubber.output.FASTQ_SCRUBBED,
    output:
        FASTA = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'output' / 'fasta.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'human_read_scrubbing' / wildcards.input_format / 'output'
    run:
        from Bio import SeqIO
        path_in = (SnakemakeUtils.load_object(Path(input.FASTQ_SCRUBBED)))[0].path
        path_out = Path(str(params.running_dir), f"{path_in.stem}.fasta")

        with path_in.open('r') as fastq_file, path_out.open('w') as fasta_file:
            # Write the FASTQ SeqRecords in FASTA format, not using SeqIO.write directly because it writes multine fasta's
            fasta_out = SeqIO.FastaIO.FastaWriter(fasta_file, wrap=None)
            fasta_out.write_file(SeqIO.parse(fastq_file, 'fastq'))
        SnakemakeUtils.dump_object([ToolIOFile(path_out)], Path(output.FASTA))


rule scrubbing_fasta_fq2fa_removed:
    """
    Converts the FASTQ containing removed reads back to FASTA.
    """
    input:
        FASTQ_REMOVED = rules.scrubbing_run_scrubber.output.FASTQ_REMOVED,
        INFORMS_SCRUBBER = rules.scrubbing_run_scrubber.output.INFORMS
    output:
        FASTA_REMOVED = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'output' / 'fasta_removed.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'human_read_scrubbing' / wildcards.input_format / 'output'
    run:
        from Bio import SeqIO

        scrubber_informs = SnakemakeUtils.load_object(Path(input.INFORMS_SCRUBBER))
        path_removed_in = (SnakemakeUtils.load_object(Path(input.FASTQ_REMOVED)))[0].path
        path_removed_out = Path(str(params.running_dir), f"{path_removed_in.stem}.fasta")

        if scrubber_informs['statistics']['count_removed'] != 0:
            with path_removed_in.open('r') as fastq_removed_file, path_removed_out.open('w') as fasta_removed_file:
                # Write the FASTQ SeqRecords in FASTA format, not using SeqIO.write directly because it writes multine fasta's
                fasta_removed_out = SeqIO.FastaIO.FastaWriter(fasta_removed_file, wrap=None)
                fasta_removed_out.write_file(SeqIO.parse(fastq_removed_file, 'fastq'))
            SnakemakeUtils.dump_object([ToolIOFile(path_removed_out)], Path(output.FASTA_REMOVED))
        else:
            SnakemakeUtils.dump_object([], Path(output.FASTA_REMOVED))


rule scrubbing_fastq_deinterleave_and_gzip:
    """
    If the input is a paired-end interleaved file, deinterleaves. Gzips in all cases.
    """
    input:
        FASTQ = Path(config['working_dir']) / human_read_scrubbing.INPUT_SCRUBBING_FASTQ,
        FASTQ_SCRUBBED = rules.scrubbing_run_scrubber.output.FASTQ_SCRUBBED,
    output:
        FASTQ_DEINTERLEAVED_GZIPPED = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'output' / 'fastq.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'human_read_scrubbing' / wildcards.input_format / 'output',
        is_paired = lambda wildcards: wildcards.input_format == 'fastq_pe'
    threads: 8
    run:
        from camel.app.components.filesystemhelper import FileSystemHelper
        from camel.app.components.files.fastqutils import FastqUtils

        # Get the FASTQ file(s)
        fastq_in = SnakemakeUtils.load_object(Path(input.FASTQ))

        if not params.is_paired:
            output.FASTQ_SINGLE_GUNZIP = input.FASTQ_SCRUBBED
            path_in = (SnakemakeUtils.load_object(Path(input.FASTQ_SCRUBBED)))[0].path
            path_out = Path(str(params.running_dir), f"{fastq_in[0].path.stem.split('.fastq')[0].split('.fq')[0]}.fastq.gz")
            FileSystemHelper.pigz_file(path_in, path_out, threads=threads)
            SnakemakeUtils.dump_object([ToolIOFile(path_out)], Path(output.FASTQ_DEINTERLEAVED_GZIPPED))
        else:
            params.running_dir.mkdir(parents=True, exist_ok=True)
            fastq_1 = Path(str(params.running_dir), f"{FastqUtils.get_sample_name(fastq_in[0].path, pattern=FastqUtils.PATTERN_FQ_SE)}_1.fastq.gz")
            fastq_2 = Path(str(params.running_dir), f"{FastqUtils.get_sample_name(fastq_in[1].path, pattern=FastqUtils.PATTERN_FQ_SE)}_2.fastq.gz")
            FastqUtils.split_interleaved_fastq((SnakemakeUtils.load_object(Path(input.FASTQ_SCRUBBED)))[0].path, fastq_1, fastq_2, gzip_output=True, pigz=True)
            SnakemakeUtils.dump_object([ToolIOFile(fastq_1), ToolIOFile(fastq_2)], Path(output.FASTQ_DEINTERLEAVED_GZIPPED))


rule scrubbing_fastq_deinterleave_and_gzip_removed:
    """
    If the input is a paired-end interleaved file, deinterleaves the removed reads. Gzips in all cases the removed reads.
    """
    input:
        FASTQ_REMOVED = rules.scrubbing_run_scrubber.output.FASTQ_REMOVED,
        INFORMS_SCRUBBER = rules.scrubbing_run_scrubber.output.INFORMS
    output:
        FASTQ_REMOVED_DEINTERLEAVED_GZIPPED = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'output' / 'fastq_removed.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'human_read_scrubbing' / wildcards.input_format / 'output',
        is_paired = lambda wildcards: wildcards.input_format == 'fastq_pe'
    run:
        from camel.app.components.filesystemhelper import FileSystemHelper
        from camel.app.components.files.fastqutils import FastqUtils

        scrubber_informs = SnakemakeUtils.load_object(Path(input.INFORMS_SCRUBBER))
        fastq_removed_in = SnakemakeUtils.load_object(Path(input.FASTQ_REMOVED))

        if scrubber_informs['statistics']['count_removed'] != 0:
            if not params.is_paired:
                path_removed_in = (SnakemakeUtils.load_object(Path(input.FASTQ_REMOVED)))[0].path
                path_removed_out = Path(str(params.running_dir), f"{fastq_removed_in[0].path.stem.split('.fastq')[0].split('.fq')[0]}.fastq.gz")
                FileSystemHelper.gzip_file(path_removed_in, path_removed_out)
                SnakemakeUtils.dump_object([ToolIOFile(path_removed_out)], Path(output.FASTQ_REMOVED_DEINTERLEAVED_GZIPPED))
            else:
                fastq_removed_1 = Path(str(params.running_dir), f"{FastqUtils.get_sample_name(fastq_removed_in[0].path, pattern=FastqUtils.PATTERN_FQ_SE)}_1.fastq.gz")
                fastq_removed_2 = Path(str(params.running_dir), f"{FastqUtils.get_sample_name(fastq_removed_in[0].path, pattern=FastqUtils.PATTERN_FQ_SE)}_2.fastq.gz")
                FastqUtils.split_interleaved_fastq((SnakemakeUtils.load_object(Path(input.FASTQ_REMOVED)))[0].path, fastq_removed_1, fastq_removed_2, gzip_output=True)
                SnakemakeUtils.dump_object([ToolIOFile(fastq_removed_1), ToolIOFile(fastq_removed_2)], Path(output.FASTQ_REMOVED_DEINTERLEAVED_GZIPPED))
        else:
            SnakemakeUtils.dump_object([],Path(output.FASTQ_REMOVED_DEINTERLEAVED_GZIPPED))

rule scrubbing_report:
    """
    Generates a tiny html report containing the tool name/version and a single phrase about how many reads were removed.
    """
    input:
        INFORMS_SCRUBBER = Path(config['working_dir']) / human_read_scrubbing.OUTPUT_SCRUBBING_INFORMS,
        REMOVED = lambda wildcards: get_removed(config, input_format=wildcards.input_format) if config['read_scrubbing'].get('export_removed_reads') else []
    output:
        HTML = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'output' / 'html.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'human_read_scrubbing' / wildcards.input_format / 'output',
        input_format = lambda wildcards: wildcards.input_format,
        export_removed_reads = config['read_scrubbing'].get('export_removed_reads')
    run:
        from camel.app.tools.ncbihumanreadscrubber.ncbihumanreadscrubberreporter import NcbiHumanReadScrubberReporter

        reporter = NcbiHumanReadScrubberReporter(Camel.get_instance())
        if params.export_removed_reads:
            SnakemakeUtils.add_pickle_inputs(reporter, input, keys=['REMOVED'])
        SnakemakeUtils.add_pickle_inputs(reporter, input, keys=['INFORMS_SCRUBBER'])
        reporter.update_parameters(input_format=str(params.input_format))
        step = Step(str(rule), reporter, Camel.get_instance(), Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule scrubbing_create_summary:
    """
    Creates the tabular summary output for the HRRT assay.
    """
    input:
        INFORMS_SCRUBBER = Path(config['working_dir']) / human_read_scrubbing.OUTPUT_SCRUBBING_INFORMS
    output:
        TSV = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'output' / 'summary_out.tsv'
    params:
        input_format = lambda wildcards: wildcards.input_format
    run:
        informs = SnakemakeUtils.load_object(Path(input.INFORMS_SCRUBBER))

        subject = 'read_pairs' if params.input_format == 'fastq_pe' else 'reads' if params.input_format == 'fastq_se' \
            else 'contigs'
        count_total = informs['statistics']['count_total']
        count_removed = informs['statistics']['count_removed']

        data_summary = {
            'scrubbing_tool_version': informs['_name'],
            f'scrubbing_{subject}_in': count_total,
            f'scrubbing_{subject}_removed': count_removed
        }

        with Path(output.TSV).open('w') as handle:
            for k, v in data_summary.items():
                handle.write('\t'.join([k, str(v)]))
                handle.write('\n')

rule scrubbing_report_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'output' / 'html-empty.io'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Human read removal', Path(output.VAL_HTML))
