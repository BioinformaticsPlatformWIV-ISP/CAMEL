import gzip
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.hybridassemblypipeline.snakefile import assembly_flye, short_read_polishing, medaka_snakemake, quality_checks

camel = Camel.get_instance()

include: assembly_flye.SNAKEFILE_FLYE
include: medaka_snakemake.SNAKEFILE_POLISHING
include: short_read_polishing.SNAKEFILE_POLISHING
include: quality_checks.SNAKEFILE_QC

#########
# Rules #
#########

rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        Path(config['working_dir'] / 'ok.txt')

rule trim_illumina:
    input:
        FQ_fwd = config['input']['illumina'][0],
        FQ_rev = config['input']['illumina'][1]
    output:
        FQ_1P = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_1P.fastq.gz",
        FQ_2P = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_2P.fastq.gz",
        FQ_1S = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_1U.fastq.gz",
        FQ_2S = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_2U.fastq.gz",
        TSV = Path(config['working_dir']) / 'trimming' / 'illumina' / 'trimming_illumina.tsv',
        HTML = Path(config['working_dir']) / 'trimming' / 'illumina' / 'fastqc_pre.html'
    params:
        dir_ = Path(config['working_dir']) / 'trimming' / 'illumina'
    threads: 4
    run:
        from camel.app.components.workflows.trimmingilluminawrapper import TrimmingIlluminaWrapper
        wrapper = TrimmingIlluminaWrapper(Path(params.dir_).absolute())
        wrapper.run_workflow([Path(input.FQ_fwd), Path(input.FQ_rev)], threads=threads)
        wrapper.output.trimmed_reads_pe[0].path.rename(Path(output.FQ_1P))
        wrapper.output.trimmed_reads_pe[1].path.rename(Path(output.FQ_2P))
        wrapper.output.trimmed_reads_se_fwd[0].path.rename(Path(output.FQ_1S))
        wrapper.output.trimmed_reads_se_rev[0].path.rename(Path(output.FQ_2S))
        wrapper.output.tsv_summary.rename(Path(output.TSV))
        wrapper.output.fastq_reports_pre[0].path.rename(Path(output.HTML))

rule trim_ont:
    input:
        FASTQ = config['input']['ont']
    output:
        FASTQ = Path(config['working_dir']) / 'trimming' / 'ont' / 'fastq.io'
    params:
        dir_ = Path(config['working_dir']) / 'trimming' / 'ont',
        filtlong_options= config.get('filtlong',{})
    threads: 4
    run:
        from camel.app.tools.filtlong.filtlong import Filtlong
        filtlong = Filtlong(camel)
        filtlong.add_input_files({'FASTQ': [ToolIOFile(Path(input.FASTQ))]})
        filtlong.update_parameters(**params.filtlong_options)
        step = Step(str(rule), filtlong, camel, params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(filtlong,output)

rule set_trimming_ont_output:
    input:
        FASTQ = rules.trim_ont.output.FASTQ
    output:
        FASTQ = Path(config['working_dir']) / 'trimming' / 'ont' / '{}_SE.fastq.gz'.format(config['name'])
    params:
        dir_ = Path(config['working_dir']) / 'trimming' / 'ont'
    run:
        input_fastq = open(SnakemakeUtils.load_object(Path(input.FASTQ))[0].path, 'rb').read()
        with gzip.open(output.FASTQ, 'wb') as handle:
            handle.write(input_fastq)

rule unicycler:
    input:
        FQ_1P = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_1P.fastq.gz",
        FQ_2P = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_2P.fastq.gz",
        FQ_SE = Path(config['working_dir']) / 'trimming' / 'ont' / '{}_SE.fastq.gz'.format(config['name'])
    output:
        FASTA = Path(config['working_dir']) / 'unicycler' / 'assembly.fasta'
    params:
        dir_ = Path(config['working_dir'])
    threads: 4
    run:
        from camel.app.tools.unicycler.unicycler import Unicycler
        unicycler_assembly = Unicycler(camel)
        unicycler_assembly.add_input_files({'FASTQ_PE':[ToolIOFile(Path(input.FQ_1P)), ToolIOFile(Path(input.FQ_2P))],
                                            'FASTQ_SE':[ToolIOFile(Path(input.FQ_SE))]})
        unicycler_assembly.update_parameters(output_dir = 'unicycler', threads=threads)
        step = Step(str(rule), unicycler_assembly, camel, params.dir_, config)
        step.run_step()
