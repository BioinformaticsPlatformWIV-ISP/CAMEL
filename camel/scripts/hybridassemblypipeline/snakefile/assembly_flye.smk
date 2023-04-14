import shutil
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step

camel = Camel.get_instance()


rule assembly_flye_run:
    """
    De-novo assembly using Flye.
    """
    input:
        FASTQ = Path(config['working_dir']) / 'trimming' / 'ont' / '{}_SE.fastq.gz'.format(config['name'])
    output:
        FASTA = Path(config['working_dir']) / 'assembly_flye' / 'flye' / 'assembly.fasta',
    params:
        running_dir = Path(config['working_dir']) / 'assembly_flye' / 'flye',
        flye_options = config.get('assembly', {}).get('flye', {}),
        read_type = 'SE'
    threads: 8
    run:
        from camel.app.tools.flye.flye import Flye
        flye = Flye(camel)
        flye.add_input_files({'FASTQ': [ToolIOFile(Path(input.FASTQ))]})
        flye.update_parameters(**params.flye_options)
        flye.update_parameters(threads=threads, output_directory=str(params.running_dir))
        step = Step(str(rule), flye, camel, params.running_dir)
        step.run_step()

rule assembly_flye_filter_contig_length:
    """
    Filters out the small contigs.
    """
    input:
        FASTA = rules.assembly_flye_run.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / 'assembly_flye' / 'filtering' / 'assembly_filtered.fasta',
    params:
        running_dir = Path(config['working_dir']) / 'assembly_flye' / 'filtering',
        min_contig_length = config['assembly'].get('min_contig_length', 0) if 'assembly' in config else 0
    run:
        from camel.app.tools.seqtk.seqtkseq import SeqtkSeq
        seqtk = SeqtkSeq(camel)
        seqtk.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA))]})
        seqtk.update_parameters(output_filename='assembly_filtered.fasta', min_length=params.min_contig_length)
        step = Step(str(rule), seqtk, camel, params.running_dir)
        step.run_step()

rule assembly_flye_quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = rules.assembly_flye_filter_contig_length.output.FASTA
    output:
        TSV = Path(config['working_dir']) / 'assembly_flye' / 'quast' / 'report.tsv'
    params:
        running_dir = Path(config['working_dir']) / 'assembly_flye' / 'quast'
    run:
        from camel.app.tools.quast.quast import Quast
        quast = Quast(camel)
        quast.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA))]})
        step = Step(str(rule), quast, camel, params.running_dir)
        step.run_step()
