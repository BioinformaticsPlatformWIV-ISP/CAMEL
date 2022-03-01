import json
import shutil
import statistics
from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import downsampling

rule downsampling_fqstats:
    """
    Determines the number of reads and bases in the input FASTQ files.
    """
    input:
        FASTQ = Path(config['working_dir']) / 'trimming_illumina' / 'input'/ 'fastq-pe.io'
    output:
        IO = Path(config['working_dir']) / 'downsampling' / 'json_fqstats.io'
    params:
        running_dir = Path(config['working_dir']) / 'downsampling'
    run:
        from camel.app.tools.fqtools.fqstats import FqStats
        fqstats = FqStats(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(fqstats, input)
        step = Step(str(rule), fqstats, Camel.get_instance(), params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_object(fqstats.informs, Path(output.IO))

rule downsampling_calculate:
    """
    Calculates the downsampling statistics.
    """
    input:
        IO = rules.downsampling_fqstats.output.IO
    output:
        JSON = Path(config['working_dir']) / 'downsampling' / 'stats.json'
    params:
        running_dir = Path(config['working_dir']) / 'downsampling',
        size_ref_genome = config['downsampling']['expected_genome_size'],
        cov_target = config['downsampling']['coverage_max']
    run:
        informs_fq = SnakemakeUtils.load_object(Path(input.IO))
        coverage_est = sum(fq['nb_of_bases'] for fq in informs_fq['stats']) / params.size_ref_genome
        downsample_factor = float(f"{params.cov_target / coverage_est:.2f}")
        data_out = {
            'total_bases': sum([fq['nb_of_bases'] for fq in informs_fq['stats']]),
            'mean_read_length': statistics.mean([
                fq['nb_of_bases'] / fq['nb_of_sequences'] for fq in informs_fq['stats']]),
            'coverage_estimated': coverage_est,
            'coverage_target': params.cov_target,
            'downsample_factor': params.cov_target / coverage_est if downsample_factor < 1 else None,
            'size_ref_genome': params.size_ref_genome,
            'nb_read_pairs': next(iter(informs_fq['stats']))['nb_of_sequences']
        }
        with open(output.JSON, 'w') as handle:
            json.dump(data_out, handle, indent=2)


rule downsampling_seqtk:
    """
    Performs downsampling with seqtk (if needed).
    """
    input:
        FASTQ_PE = Path(config['working_dir']) / 'trimming_illumina' / 'input' / 'fastq-pe.io',
        JSON = rules.downsampling_calculate.output.JSON
    output:
        FASTQ_PE = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_FASTQ_PE,
        INFORMS = Path(config['working_dir']) / 'downsampling' / 'informs_seqtk.io'
    params:
        dir_working = Path(config['working_dir']) / 'downsampling'
    run:
        import logging
        from camel.app.tools.seqtk.seqtksubsample import SeqtkSubsample

        # Parse statistics
        with open(input.JSON) as handle:
            data_ds = json.load(handle)
        if data_ds['downsample_factor'] is None:
            logging.info("No downsampling required, skipping seqtk")
            shutil.copyfile(input.FASTQ_PE, output.FASTQ_PE)
            SnakemakeUtils.dump_object(None, Path(output.INFORMS))
        else:
            seqtk = SeqtkSubsample(Camel.get_instance())
            seqtk.add_input_files({'FASTQ_PE': SnakemakeUtils.load_object(Path(input.FASTQ_PE))})
            step = Step(str(rule), seqtk, Camel.get_instance(), Path(params.dir_working))
            seqtk.update_parameters(fraction=float(f"{data_ds['downsample_factor']:.2f}"))
            step.run_step()
            SnakemakeUtils.dump_object(seqtk.tool_outputs['FASTQ_PE'], Path(output.FASTQ_PE))
            SnakemakeUtils.dump_object(seqtk.informs, Path(output.INFORMS))

rule downsampling_report:
    """
    Creates the downsampling report.
    """
    input:
        JSON = rules.downsampling_calculate.output.JSON,
        INFORMS_seqtk = rules.downsampling_seqtk.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_REPORT
    params:
        dir_working = Path(config['working_dir']) / 'downsampling'
    run:
        from camel.app.tools.pipelines.downsampling.reporterdownsampling import ReporterDownsampling
        reporter = ReporterDownsampling(Camel.get_instance())
        with open(input.JSON) as handle:
            stats = json.load(handle)
        reporter.add_input_informs({'stats': stats})
        informs_seqtk = SnakemakeUtils.load_object(Path(input.INFORMS_seqtk))
        if informs_seqtk is not None:
            reporter.add_input_informs({'seqtk': informs_seqtk})
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_working)
        step.run_step()
        SnakemakeUtils.dump_object(reporter.tool_outputs['HTML'], Path(output.HTML))
