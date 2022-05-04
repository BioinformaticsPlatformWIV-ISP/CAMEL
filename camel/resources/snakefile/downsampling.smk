import json
import shutil
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
        FASTQ = Path(config['working_dir']) / downsampling.INPUT_DOWNSAMPLING_FASTQ
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
        cov_target = config['downsampling']['coverage_max'],
        is_paired = config.get('read_type', 'illumina') == 'illumina'
    run:
        import statistics
        informs_fq = SnakemakeUtils.load_object(Path(input.IO))
        coverage_est = sum(fq['nb_of_bases'] for fq in informs_fq['stats']) / params.size_ref_genome
        downsample_factor = float(f"{params.cov_target / coverage_est:.2f}")
        key_nb_reads = 'nb_read_pairs' if params.is_paired else 'nb_reads'
        data_out = {
            'total_bases': sum([fq['nb_of_bases'] for fq in informs_fq['stats']]),
            'mean_read_length': statistics.mean([
                fq['nb_of_bases'] / fq['nb_of_sequences'] for fq in informs_fq['stats']]),
            'coverage_estimated': coverage_est,
            'coverage_target': params.cov_target,
            'downsample_factor': params.cov_target / coverage_est if downsample_factor < 1 else None,
            'size_ref_genome': params.size_ref_genome,
            f'{key_nb_reads}_in': next(iter(informs_fq['stats']))['nb_of_sequences']
        }
        with open(output.JSON, 'w') as handle:
            json.dump(data_out, handle, indent=2)


rule downsampling_seqtk:
    """
    Performs downsampling with seqtk (if needed).
    """
    input:
        FASTQ = Path(config['working_dir']) / downsampling.INPUT_DOWNSAMPLING_FASTQ,
        JSON = rules.downsampling_calculate.output.JSON
    output:
        FASTQ = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_FASTQ,
        INFORMS = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_INFORMS
    params:
        dir_working = Path(config['working_dir']) / 'downsampling',
        is_paired = config.get('read_type','illumina') == 'illumina'
    run:
        import logging
        from camel.app.tools.seqtk.seqtksubsample import SeqtkSubsample

        # Parse statistics
        with open(input.JSON) as handle:
            data_ds = json.load(handle)
        if data_ds['downsample_factor'] is None:
            logging.info("No downsampling required, skipping seqtk")
            shutil.copyfile(input.FASTQ, output.FASTQ)
            SnakemakeUtils.dump_object([], Path(output.INFORMS))
        else:
            seqtk = SeqtkSubsample(Camel.get_instance())
            key_fastq = 'FASTQ_PE' if params.is_paired else 'FASTQ'
            seqtk.add_input_files({key_fastq: SnakemakeUtils.load_object(Path(input.FASTQ))})
            step = Step(str(rule), seqtk, Camel.get_instance(), Path(params.dir_working))
            seqtk.update_parameters(fraction=float(f"{data_ds['downsample_factor']:.2f}"))
            step.run_step()
            SnakemakeUtils.dump_object(seqtk.tool_outputs[key_fastq], Path(output.FASTQ))
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
        dir_working = Path(config['working_dir']) / 'downsampling',
        is_paired = config.get('read_type','illumina') == 'illumina'
    run:
        from camel.app.tools.pipelines.downsampling.reporterdownsampling import ReporterDownsampling
        reporter = ReporterDownsampling(Camel.get_instance())
        with open(input.JSON) as handle:
            stats = json.load(handle)
        reporter.add_input_informs({'stats': stats})
        if params.is_paired:
            reporter.update_parameters(is_paired=None)
        else:
            reporter.update_parameters(is_paired=False)
        informs_seqtk = SnakemakeUtils.load_object(Path(input.INFORMS_seqtk))
        if len(informs_seqtk) > 0:
            reporter.add_input_informs({'seqtk': informs_seqtk})
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_working)
        step.run_step()
        SnakemakeUtils.dump_object(reporter.tool_outputs['HTML'], Path(output.HTML))


rule downsampling_summary_out:
    """
    Exports the summary information for the downsampling workflow.
    """
    input:
        JSON = rules.downsampling_calculate.output.JSON,
        INFORMS_seqtk = rules.downsampling_seqtk.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_SUMMARY
    params:
        is_paired = config.get('read_type', 'illumina') == 'illumina'
    run:
        # Parse calculate output
        with open(input.JSON) as handle:
            data_ds = json.load(handle)

        # Add the number of output reads / read pairs to the summary output
        key_nb_reads = 'nb_read_pairs' if params.is_paired else 'nb_reads'
        key_out = f'{key_nb_reads}_out'
        if data_ds['downsample_factor'] is None:
            data_ds[key_out] = data_ds[f'{key_nb_reads}_in']
        else:
            informs_in = SnakemakeUtils.load_object(Path(input.INFORMS_seqtk))
            if params.is_paired:
                data_ds[key_out] = informs_in['reads_count'] // 2
            else:
                data_ds[key_out] = informs_in['reads_count']

        # Create output
        with open(output.TSV, 'w') as handle:
            for k, v in sorted(data_ds.items()):
                handle.write('\t'.join([
                    f'downsampling_{k}', f'{v:.2f}' if isinstance(v, float) else str(v)
                ]))
                handle.write('\n')
