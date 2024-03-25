from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.pipelineexecutionerror import PipelineExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import human_read_scrubbing

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
    run:
        from camel.app.components.filesystemhelper import FileSystemHelper
        from camel.app.components.files.fastqutils import FastqUtils

        # Get the FASTQ file(s)
        fastq_in = SnakemakeUtils.load_object(Path(input.FASTQ))
        nb_of_fq_files = len(fastq_in)

        if nb_of_fq_files == 1:
            path_in = fastq_in[0].path
            path_out = Path(str(params.running_dir), path_in.name.replace('.gz', ''))
            FileSystemHelper.gzip_extract(path_in, path_out)
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
        INFORMS = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'scrubbing' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'human_read_scrubbing' / wildcards.input_format / 'scrubbing',
        input_format = lambda wildcards: wildcards.input_format
    run:
        from camel.app.tools.ncbihumanreadscrubber.ncbihumanreadscrubber import NcbiHumanReadScrubber
        if params.input_format != 'fasta':
            # Get the FASTQ file(s)
            fastq_in = SnakemakeUtils.load_object(Path(input.FASTQ))
            fqfile_number = len(fastq_in)
            interleaved = 'true' if fqfile_number == 2 else 'false'
        else:
            interleaved = 'false'
        scrubber = NcbiHumanReadScrubber(camel)
        step = Step(str(rule), scrubber, camel, Path(str(params.running_dir)))
        outputfile_scrubbing = str(Path(output.FASTQ_SCRUBBED).with_suffix('.fastq')) if params.input_format != 'fasta' \
            else str(Path(str(params.running_dir), (SnakemakeUtils.load_object(Path(input.FASTQ_SINGLE_GUNZIP)))[0].path.name))
        scrubber.update_parameters(interleaved=interleaved, outputfile=outputfile_scrubbing)
        SnakemakeUtils.add_pickle_inputs(scrubber, input, excluded_keys=['FASTQ'])
        step.run_step()
        if config['input_type'] == 'hybrid':
            if params.input_format == 'fastq_pe':
                scrubber.informs['_tag'] = 'Illumina'
            else:
                scrubber.informs['_tag'] = 'ONT'
        SnakemakeUtils.dump_tool_outputs(scrubber, output, keys=['FASTQ_SCRUBBED', 'INFORMS'])

rule scrubbing_fasta_fq2fa:
    """
    Convert the FASTQ back to FASTA in order for the rest of the pipeline to be able 
    to use it as output or as input in bacterial pipelines.
    """
    input:
        FASTQ_SCRUBBED = rules.scrubbing_run_scrubber.output.FASTQ_SCRUBBED
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

rule scrubbing_fastq_deinterleave_and_gzip:
    """
    If the input is a paired-end interleaved file, deinterleaves. Gzips in all cases.
    """
    input:
        FASTQ = Path(config['working_dir']) / human_read_scrubbing.INPUT_SCRUBBING_FASTQ,
        FASTQ_SCRUBBED = rules.scrubbing_run_scrubber.output.FASTQ_SCRUBBED
    output:
        FASTQ_DEINTERLEAVED_GZIPPED = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'output' / 'fastq.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'human_read_scrubbing' / wildcards.input_format / 'output',
        is_paired = lambda wildcards: wildcards.input_format == 'fastq_pe'
    run:
        from camel.app.components.filesystemhelper import FileSystemHelper
        from camel.app.components.files.fastqutils import FastqUtils

        # Get the FASTQ file(s)
        fastq_in = SnakemakeUtils.load_object(Path(input.FASTQ))
        if not params.is_paired:
            output.FASTQ_SINGLE_GUNZIP = input.FASTQ_SCRUBBED
            path_in = (SnakemakeUtils.load_object(Path(input.FASTQ_SCRUBBED)))[0].path
            path_out = Path(str(params.running_dir), f"{fastq_in[0].path.stem.split('.fastq')[0].split('.fq')[0]}.fastq.gz")
            FileSystemHelper.gzip_file(path_in, path_out)
            SnakemakeUtils.dump_object([ToolIOFile(path_out)], Path(output.FASTQ_DEINTERLEAVED_GZIPPED))
        else:
            params.running_dir.mkdir(parents=True, exist_ok=True)
            fastq_1 = Path(str(params.running_dir), f"{FastqUtils.get_sample_name(fastq_in[0].path, pattern=FastqUtils.PATTERN_FQ_SE)}_1.fastq.gz")
            fastq_2 = Path(str(params.running_dir), f"{FastqUtils.get_sample_name(fastq_in[1].path, pattern=FastqUtils.PATTERN_FQ_SE)}_2.fastq.gz")
            FastqUtils.split_interleaved_fastq((SnakemakeUtils.load_object(Path(input.FASTQ_SCRUBBED)))[0].path, fastq_1,  fastq_2, gzip_output=True)
            SnakemakeUtils.dump_object([ToolIOFile(fastq_1), ToolIOFile(fastq_2)], Path(output.FASTQ_DEINTERLEAVED_GZIPPED))

rule scrubbing_report:
    """
    Generates a tiny html report containing the tool name/version and a single phrase about how many reads were removed.
    """
    input:
        INFORMS_tools = Path(config['working_dir']) / human_read_scrubbing.OUTPUT_SCRUBBING_INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'output' / 'html.io',
        TSV = Path(config['working_dir']) / 'human_read_scrubbing' / '{input_format}' / 'output' / 'summary_out.tsv'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'human_read_scrubbing' / wildcards.input_format / 'output',
        input_format = lambda wildcards: wildcards.input_format
    run:
        from camel.app.components.html.htmlreportsection import HtmlReportSection

        # Parse informs
        hrrt_informs = SnakemakeUtils.load_object(Path(input.INFORMS_tools))
        count_removed = hrrt_informs['statistics']['count_removed']
        count_in = hrrt_informs['statistics']['count_total']
        if hrrt_informs['statistics']['count_removed'] == hrrt_informs['statistics']['count_total']:
            raise PipelineExecutionError(
                'All reads/contigs were removed from the input file(s) during scrubbing. If this is not expected, '
                'try disabling the human read scrubbing step.')

        # Create the report section
        section = HtmlReportSection('Human read removal', subtitle=hrrt_informs['_name'])
        subject = 'read_pairs' if params.input_format == 'fastq_pe' else 'reads' if params.input_format == 'fastq_se' else 'contigs'
        section.add_table([
            [f'Total {subject.replace("_", " ")}', f'{count_in:,}'],
            [f'Removed {subject.replace("_", " ")}', f'{count_removed:,}'],
            [f'Removed %', f'{100 * count_removed / count_in:.2f}'],
        ],['Category', 'Number'],[('class', 'data')])
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.VAL_HTML))

        # Create the summary output
        data_summary = {'scrubbing_tool_version': hrrt_informs['_name'], f'scrubbing_{subject}_in': count_in,
                        f'scrubbing_{subject}_out': count_removed}
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
