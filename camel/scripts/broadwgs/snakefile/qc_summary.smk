from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.broadwgs.snakefile import alignment, bam_to_cram, variant_calling

camel = Camel.get_instance()


rule generate_qc_summary:
    input:
        mark_duplicates_metrics = Path(config['working_dir']) / 'qc' / 'mark_duplicates' / "duplicate_metrics.txt.io",
        alignment_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f'{config["sample"]}.agg.alignment_summary_metrics',
        insert_size_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f'{config["sample"]}.agg.insert_size_metrics',
        wgs_metrics = Path(config['working_dir']) / "qc" / "wgs_metrics" / f'{config["sample"]}.wgs.metrics.txt',
        variant_calling_metrics = Path(config['working_dir']) / "qc" / "variant_calling_metrics" / f'{config["sample"]}.variant_calling_metrics.io'
    output:
        QC_summary = Path(config['working_dir']) / "qc" / "QC_summary.txt"
    run:
        alignment_metrics_threshold = {
            'percent_duplication' : {'val' : 0.1, 'threshold' : 'upper'},
            'percent_chimeras' : {'val': 0.01, 'threshold' : 'upper'},
            'median_insert_size' : {'val': 300, 'threshold': 'lower'},
            'mean_coverage' : {'val': 30, 'threshold': 'lower'},
            'evenness_coverage' : {'val': 1, 'threshold' : 'both'}
        }

        variant_calling_metrics_threshold = {
            "TITV_DBSNP" : {'val': 2, 'threshold': 'both'},
            "TITV_NOVEL" : {'val': 2, 'threshold': 'both'}
        }

        alignment_metrics = {}
        variant_calling_metrics = {}

        # Mark Duplicates metrics
        mark_duplicates_file = Path(SnakemakeUtils.load_object(Path(input.mark_duplicates_metrics))[0].path)
        with open(mark_duplicates_file, 'r') as handle:
            lines = handle.readlines()
            metrics_class_index = lines.index("## METRICS CLASS\tpicard.sam.DuplicationMetrics\n")

        header = lines[metrics_class_index + 1].split("\t")
        values = lines[metrics_class_index + 2].split("\t")

        mark_duplicates = dict(zip(header, values))
        alignment_metrics['percent_duplication'] = mark_duplicates['PERCENT_DUPLICATION']

        # Alignment summary metrics
        with open(input.alignment_summary_metrics, 'r') as handle:
            lines = handle.readlines()
            metrics_class_index = lines.index("## METRICS CLASS\tpicard.analysis.AlignmentSummaryMetrics\n")

        header = lines[metrics_class_index + 1].split("\t")
        values = lines[metrics_class_index + 4].split("\t")

        agg_alignment_summary = dict(zip(header, values))
        alignment_metrics['percent_chimeras'] = agg_alignment_summary['PCT_CHIMERAS']

        # Insert size metrics
        with open(input.insert_size_metrics, 'r') as handle:
            lines = handle.readlines()
            metrics_class_index = lines.index("## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\n")

        header = lines[metrics_class_index + 1].split("\t")
        values = lines[metrics_class_index + 2].split("\t")

        agg_insert_size = dict(zip(header, values))
        alignment_metrics['median_insert_size'] = agg_insert_size['MEDIAN_INSERT_SIZE']

        # WGS metrics
        with open(input.wgs_metrics, 'r') as handle:
            lines = handle.readlines()
            metrics_class_index = lines.index("## METRICS CLASS\tpicard.analysis.WgsMetrics\n")
            histogram_start_index = lines.index("## HISTOGRAM\tjava.lang.Integer\n")

        header = lines[metrics_class_index + 1].split("\t")
        values = lines[metrics_class_index + 2].split("\t")

        wgs_metrics = dict(zip(header, values))
        mean_coverage = float(wgs_metrics['MEAN_COVERAGE'])
        genome_territory = int(wgs_metrics['GENOME_TERRITORY'])
        lowest_20 = float(int(genome_territory) * 0.2)
        #
        # coverage_histo = lines[histogram_start_index+3:-1]
        # coverage_histo = [l.rstrip().split("\t") for l in coverage_histo]
        #
        # cov_values = [int(l[1]) for l in coverage_histo]
        # histo_cumsum = [sum(cov_values[:i+1]) for i, value in enumerate(cov_values)]
        #
        # tmp = histo_cumsum + [int(lowest_20)]
        # tmp.sort()
        # find_index = tmp.index(int(lowest_20))
        # val_before = histo_cumsum[find_index - 1]
        # val_after = histo_cumsum[find_index]
        #
        # cov_20 = (((val_after - lowest_20) / (val_after - val_before)) + find_index - 1)
        #
        alignment_metrics['mean_coverage'] = mean_coverage
        # fold80 = mean_coverage / cov_20
        # alignment_metrics['evenness_coverage'] = fold80

        # Variant calling metrics
        var_calling_file = Path(SnakemakeUtils.load_object(Path(input.variant_calling_metrics))[0].path)
        with open(var_calling_file, 'r') as handle:
            lines = handle.readlines()
            metrics_class_index = lines.index("## METRICS CLASS\tpicard.vcf.CollectVariantCallingMetrics$VariantCallingSummaryMetrics\n")

        header = lines[metrics_class_index + 1].split("\t")
        values = lines[metrics_class_index + 2].split("\t")

        variant_calling_summary = dict(zip(header, values))
        variant_calling_metrics["TITV_DBSNP"] = variant_calling_summary['DBSNP_TITV']
        variant_calling_metrics["TITV_NOVEL"] = variant_calling_summary['NOVEL_TITV']

        # Print all to output file
        alignment_metrics_metrics_names = ['mean_coverage', 'percent_duplication', 'median_insert_size', 'percent_chimeras']#'evenness_coverage',
        variant_calling_metrics_names = ['TITV_DBSNP', 'TITV_NOVEL']
        with open(output.QC_summary, 'w') as outhandle:
            for aln_metric in alignment_metrics_metrics_names:
                outhandle.write("\t".join([aln_metric, str(alignment_metrics[aln_metric]),
                                           str(alignment_metrics_threshold[aln_metric]['val']),
                                           alignment_metrics_threshold[aln_metric]['threshold']]) + "\n")

            for var_metric in variant_calling_metrics_names:
                outhandle.write("\t".join([var_metric, str(variant_calling_metrics[var_metric]),
                                           str(variant_calling_metrics_threshold[var_metric]['val']),
                                           variant_calling_metrics_threshold[var_metric]['threshold']]) + "\n")
