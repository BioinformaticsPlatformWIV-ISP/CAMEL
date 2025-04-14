from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import downsampling


rule downsampling_read_stats:
    """
    Determines the number of reads and bases in the input FASTQ files.
    """
    input:
        FASTQ = Path(config['working_dir']) / downsampling.INPUT_DOWNSAMPLING_FASTQ
    output:
        INFORMS = Path(config['working_dir']) / 'downsampling' / '{read_key}' / 'informs-fqstats.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'downsampling' / wildcards.read_key
    run:
        from camel.app.tools.seqtk.seqtksize import SeqtkSize
        seqtk_size = SeqtkSize(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(seqtk_size, input)
        step = Step(str(rule), seqtk_size, Camel.get_instance(), Path(str(params.dir_)))
        step.run_step()
        SnakemakeUtils.dump_object(seqtk_size.informs, Path(output.INFORMS))

rule downsampling_calculate:
    """
    Calculates if downsampling is required and associated statistics.
    """
    input:
        INFORMS_stats = rules.downsampling_read_stats.output.INFORMS
    output:
        JSON = Path(config['working_dir']) / 'downsampling' / '{read_key}' / 'calculate_stats' / 'json.io',
        INFORMS = Path(config['working_dir']) / 'downsampling' / '{read_key}' / 'calculate_stats' / 'informs.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'downsampling' / wildcards.read_key / 'calculate_stats',
        fasta_ref = config['reference'].get('fasta'),
        cov_target = config['downsampling']['coverage_max'],
        is_paired = lambda wildcards: wildcards.read_key == 'fastq_pe'
    run:
        from camel.app.components.files.fastautils import FastaUtils
        from camel.app.tools.pipelines.downsampling.downsamplecalculation import DownsampleCalculation
        ds_calc = DownsampleCalculation(Camel.get_instance())
        step = Step(str(rule), ds_calc, Camel.get_instance(), Path(str(params.dir_)))

        if params.fasta_ref is None:
            raise KeyError('Reference FASTA file is required to estimate genome size')
        size_ref_genome = FastaUtils.count_bases(Path(params.fasta_ref))

        SnakemakeUtils.add_pickle_inputs(ds_calc, input)
        ds_calc.update_parameters(size_ref_genome=size_ref_genome, cov_target=params.cov_target)
        if params.is_paired is True:
            ds_calc.update_parameters(is_paired=None)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(ds_calc, output)

rule downsampling_seqtk:
    """
    Performs downsampling with seqtk (if needed).
    """
    input:
        FASTQ = Path(config['working_dir']) / downsampling.INPUT_DOWNSAMPLING_FASTQ,
        INFORMS = rules.downsampling_calculate.output.INFORMS
    output:
        FASTQ = Path(config['working_dir']) / 'downsampling' / '{read_key}' / 'seqtk' / 'fastq.io',
        INFORMS = Path(config['working_dir']) / 'downsampling' / '{read_key}' / 'seqtk' / 'informs.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'downsampling' / wildcards.read_key / 'seqtk',
        is_paired = lambda wildcards: wildcards.read_key == 'fastq_pe'
    run:
        import logging
        import shutil
        from camel.app.tools.seqtk.seqtksubsample import SeqtkSubsample

        # Create working directory
        dir_ = Path(str(params.dir_)).absolute()
        dir_.mkdir(exist_ok=True, parents=True)

        # Parse statistics
        data_ds = SnakemakeUtils.load_object(Path(input.INFORMS))['stats']
        if data_ds['downsample_factor'] is None:
            # No downsampling -> create empty output files
            logging.info("No downsampling required, skipping seqtk")
            shutil.copyfile(input.FASTQ, output.FASTQ)
            SnakemakeUtils.dump_object([], Path(output.INFORMS))
        else:
            # Downsampling -> run seqtk
            seqtk = SeqtkSubsample(Camel.get_instance())
            key_fastq = 'FASTQ_PE' if params.is_paired else 'FASTQ'
            seqtk.add_input_files({key_fastq: SnakemakeUtils.load_object(Path(input.FASTQ))})
            step = Step(str(rule), seqtk, Camel.get_instance(), dir_)
            seqtk.update_parameters(fraction=float(f"{data_ds['downsample_factor']:.6f}"))
            step.run_step()
            SnakemakeUtils.dump_object(seqtk.tool_outputs[key_fastq], Path(output.FASTQ))
            SnakemakeUtils.dump_object(seqtk.informs, Path(output.INFORMS))

rule downsampling_report:
    """
    Creates the downsampling report.
    """
    input:
        INFORMS_stats = rules.downsampling_calculate.output.INFORMS,
        INFORMS_seqtk = rules.downsampling_seqtk.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / 'downsampling' / '{read_key}' / 'html.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'downsampling' / wildcards.read_key,
        is_paired = lambda wildcards: wildcards.read_key == 'fastq_pe'
    run:
        from camel.app.tools.pipelines.downsampling.reporterdownsampling import ReporterDownsampling
        reporter = ReporterDownsampling(Camel.get_instance())
        reporter.add_input_informs({'stats': SnakemakeUtils.load_object(Path(input.INFORMS_stats))})
        if params.is_paired:
            reporter.update_parameters(is_paired=None)
        else:
            reporter.update_parameters(is_paired=False)
        informs_seqtk = SnakemakeUtils.load_object(Path(input.INFORMS_seqtk))
        if len(informs_seqtk) > 0:
            reporter.add_input_informs({'seqtk': informs_seqtk})
        step = Step(str(rule), reporter, Camel.get_instance(), Path(str(params.dir_)))
        step.run_step()
        SnakemakeUtils.dump_object(reporter.tool_outputs['HTML'], Path(output.HTML))

rule downsampling_summary_out:
    """
    Exports the summary information for the downsampling workflow.
    """
    input:
        INFORMS_stats = rules.downsampling_calculate.output.INFORMS,
        INFORMS_seqtk = rules.downsampling_seqtk.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / 'downsampling' / '{read_key}' / 'summary_downsampling.tsv'
    params:
        is_paired = lambda wildcards: wildcards.read_key == 'fastq_pe',
        read_key = lambda wildcards: wildcards.read_key,
        add_read_key_prefix = config['input_type'] == 'hybrid'
    run:
        # Parse calculation output
        data_ds = SnakemakeUtils.load_object(Path(input.INFORMS_stats))['stats']

        # Add the number of output reads / read pairs to the summary output
        key_nb_reads = 'nb_read_pairs' if params.is_paired else 'nb_reads'
        key_out = f'{key_nb_reads}_out'
        if data_ds['downsample_factor'] is None:
            data_ds[key_out] = data_ds[f'{key_nb_reads}_in']
            data_ds['downsample_factor'] = '-'
        else:
            informs_in = SnakemakeUtils.load_object(Path(input.INFORMS_seqtk))
            if params.is_paired:
                data_ds[key_out] = informs_in['reads_count'] // 2
            else:
                data_ds[key_out] = informs_in['reads_count']

        # Create output
        with open(output.TSV, 'w') as handle:
            for k, v in sorted(data_ds.items()):
                # Add read key information when multiple FASTQ datasets were processed
                key = f'downsampling_{k}' if not params.add_read_key_prefix else f'downsampling_{params.read_key}_{k}'
                handle.write('\t'.join([key, f'{v:.2f}' if isinstance(v, float) else str(v)]))
                handle.write('\n')
