from pathlib import Path

from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.files.fastqutils import FastqUtils
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
        FASTA = Path(config['working_dir']) / 'human_read_scrubbing' / 'input' / 'fasta.io'
    output:
        FASTQ_from_fasta = Path(config['working_dir']) / 'human_read_scrubbing' / 'input' / 'fasta2fastq.io'
    params:
        running_dir = Path(config['working_dir']) / 'human_read_scrubbing' / 'input'
    run:
        from camel.app.components.files.fastautils import FastaUtils

        fasta_path_in = (SnakemakeUtils.load_object(Path(input.FASTA)))[0].path
        fastq_path_out = params.running_dir / f"{fasta_path_in.stem}.fastq"

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
        FASTQ_SINGLE_GUNZIP = Path(config['working_dir']) / 'human_read_scrubbing' / 'input' / 'fastq_gunzip_interleaved.io'
    params:
        running_dir = Path(config['working_dir']) / 'human_read_scrubbing' / 'input'
    run:
        # Get the FASTQ file(s)
        fastq_in = SnakemakeUtils.load_object(Path(input.FASTQ))
        fqfile_number = len(fastq_in)

        if fqfile_number == 1:
            path_in = fastq_in[0].path
            path_out = params.running_dir / path_in.name.replace('.gz', '')
            FileSystemHelper.gzip_extract(path_in, path_out)
            SnakemakeUtils.dump_object([ToolIOFile(path_out)], Path(output.FASTQ_SINGLE_GUNZIP))
        else:
            interleaved_out = params.running_dir / f"{fastq_in[0].path.stem}_interleaved.fastq"
            FastqUtils.convert_fastqs_to_interleaved_fastq(fastq_in[0].path, fastq_in[1].path, params.running_dir / f"{fastq_in[0].path.stem}_interleaved.fastq")
            SnakemakeUtils.dump_object([ToolIOFile(interleaved_out)], Path(output.FASTQ_SINGLE_GUNZIP))


rule scrubbing_run_scrubber:
    """
    Runs the NCBI human read scrubber on the input FASTQ, only accepts single gunzipped FASTQ files.
    """
    input:
        FASTQ = Path(config['working_dir']) / human_read_scrubbing.INPUT_SCRUBBING_FASTQ if not 'fasta' in config['input'] else [],
        FASTQ_SINGLE_GUNZIP = rules.scrubbing_fastq_interleave_and_gunzip.output.FASTQ_SINGLE_GUNZIP if not 'fasta' in config['input'] else rules.scrubbing_fasta_fa2fq.output.FASTQ_from_fasta
    output:
        FASTQ_SCRUBBED = Path(config['working_dir']) / 'human_read_scrubbing' / 'scrubbing' / 'fastq_scrubbed.io',
        INFORMS = Path(config['working_dir']) / human_read_scrubbing.OUTPUT_SCRUBBING_INFORMS,
    params:
        running_dir = Path(config['working_dir']) / 'human_read_scrubbing' / 'scrubbing'
    run:
        from camel.app.tools.ncbihumanreadscrubber.ncbihumanreadscrubber import NcbiHumanReadScrubber
        if not 'fasta' in config['input']:
            # Get the FASTQ file(s)
            fastq_in = SnakemakeUtils.load_object(Path(input.FASTQ))
            fqfile_number = len(fastq_in)
            interleaved = 'true' if fqfile_number == 2 else 'false'
        else:
            interleaved = 'false'
        scrubber = NcbiHumanReadScrubber(camel)
        step = Step(str(rule), scrubber, camel, params.running_dir)
        scrubber.update_parameters(interleaved=interleaved, outputfile=(str(Path(output.FASTQ_SCRUBBED).with_suffix('.fastq')) if not 'fasta' in config['input'] else
                                                                        str(params.running_dir / (SnakemakeUtils.load_object(Path(input.FASTQ_SINGLE_GUNZIP)))[0].path.name)))
        SnakemakeUtils.add_pickle_inputs(scrubber, input, excluded_keys=['FASTQ'])
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(scrubber, output, keys=['FASTQ_SCRUBBED', 'INFORMS'])


rule scrubbing_fasta_fq2fa:
    """
    Convert the FASTQ back to FASTA in order for the rest of the pipeline to be able 
    to use it as output or as input in bacterial pipelines.
    """
    input:
        FASTQ_SCRUBBED = rules.scrubbing_run_scrubber.output.FASTQ_SCRUBBED
    output:
        FASTA = Path(config['working_dir']) / human_read_scrubbing.OUTPUT_SCRUBBING_FASTA
    params:
        running_dir = Path(config['working_dir']) / 'human_read_scrubbing' / 'output'
    run:
        path_in = (SnakemakeUtils.load_object(Path(input.FASTQ_SCRUBBED)))[0].path
        path_out = params.running_dir / f"{path_in.stem}.fasta"

        with path_in.open('r') as fastq_file, path_out.open('w') as fasta_file:
            for record in SeqIO.parse(fastq_file, 'fastq'):
                # Write the SeqRecord in FASTA format
                fasta_file.write(f">{record.id}\n{record.seq}\n")
        SnakemakeUtils.dump_object([ToolIOFile(path_out)], Path(output.FASTA))


rule scrubbing_fastq_deinterleave_and_gzip:
    """
    If the input is a paired-end interleaved file, deinterleaves. Gzips in all cases.
    """
    input:
        FASTQ = Path(config['working_dir']) / human_read_scrubbing.INPUT_SCRUBBING_FASTQ,
        FASTQ_SCRUBBED = rules.scrubbing_run_scrubber.output.FASTQ_SCRUBBED
    output:
        FASTQ_DEINTERLEAVED_GZIPPED = Path(config['working_dir']) / human_read_scrubbing.OUTPUT_SCRUBBING_FASTQ
    params:
        running_dir = Path(config['working_dir']) / 'human_read_scrubbing' / 'output'
    run:
        # Get the FASTQ file(s)
        fastq_in = SnakemakeUtils.load_object(Path(input.FASTQ))
        fqfile_number = len(fastq_in)
        if fqfile_number == 1:
            output.FASTQ_SINGLE_GUNZIP = input.FASTQ_SCRUBBED
            path_in = (SnakemakeUtils.load_object(Path(input.FASTQ_SCRUBBED)))[0].path
            path_out = params.running_dir / f"{fastq_in[0].path.stem.split('.fastq')[0].split('.fq')[0]}.fastq.gz"
            command = Command(f'gzip -c {path_in} > {path_out}')
            command.run(path_out.parent)
            if not command.returncode == 0:
                raise PipelineExecutionError(f"Cannot unzip input file: {command.stderr}")
            SnakemakeUtils.dump_object([ToolIOFile(path_out)], Path(output.FASTQ_DEINTERLEAVED_GZIPPED))
        else:
            params.running_dir.mkdir(parents=True, exist_ok=True)
            fastq_1 = params.running_dir / f"{FastqUtils.get_sample_name(fastq_in[0].path, pattern=FastqUtils.PATTERN_FQ_SE)}.fastq.gz"
            fastq_2 = params.running_dir / f"{FastqUtils.get_sample_name(fastq_in[1].path, pattern=FastqUtils.PATTERN_FQ_SE)}.fastq.gz"
            FastqUtils.split_interleaved_fastq((SnakemakeUtils.load_object(Path(input.FASTQ_SCRUBBED)))[0].path, fastq_1,  fastq_2, gzip_output=True)
            SnakemakeUtils.dump_object([ToolIOFile(fastq_1), ToolIOFile(fastq_2)], Path(output.FASTQ_DEINTERLEAVED_GZIPPED))


rule scrubbing_report:
    """
    Generates a tiny html report containing the tool name/version and a single phrase about how many reads were removed.
    """
    input:
        INFORMS_tools = Path(config['working_dir']) / human_read_scrubbing.OUTPUT_SCRUBBING_INFORMS
    output:
        # JSON = Path(config['working_dir']) / human_read_scrubbing.OUTPUT_SCRUBBING_SUMMARY_JSON, # todo json output for hera
        VAL_HTML = Path(config['working_dir']) / human_read_scrubbing.OUTPUT_SCRUBBING_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'human_read_scrubbing' / 'output'
    run:
        from camel.app.components.html.htmlreportsection import HtmlReportSection

        hrrt_informs = SnakemakeUtils.load_object(Path(input.INFORMS_tools))
        count_removed = hrrt_informs['statistics']['count_removed']
        count_total = hrrt_informs['statistics']['count_total']
        if count_removed == count_total:
            raise PipelineExecutionError('All reads/contigs were removed from the input file(s) during scrubbing. If this is not expected, try disabling the human read scrubbing step.')
        section = HtmlReportSection('Human Read Removal', subtitle=hrrt_informs['_name'])
        section.add_paragraph(f'Removed {count_removed:,} out of {count_total:,} reads/contigs.')
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.VAL_HTML))
