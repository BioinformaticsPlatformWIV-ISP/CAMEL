# Functions to parse some of the metrics files from the WGS pipeline to generate the QC summary file

def parse_metric_file(metric_file, header_str, metric):
    """
    INPUT:
    - metric_file: path to the metric file
    - header_str:
        Picard metrics files contain the following information, preceeded by a header line:
        - Command
        - Date/time stamp
        - Summary Table containing the metrics
        - Histogram
        We need information from the summary table. To get the index of the first line of this table, we look forthe index
        of the specific header line (header_str).
    - metric: Exact metric we want to read from the summary table

    OUTPUT:
    Returns the metric
    """
    with open(metric_file, 'r') as handle:
        lines = handle.readlines()
        metrics_class_index = lines.index(header_str)

    header = lines[metrics_class_index + 1].split("\t")
    values = lines[metrics_class_index + 2].split("\t")

    parsed_metrics = dict(zip(header, values))

    return parsed_metrics[metric]

def get_fold_80(coverage_file, mean_coverage, percentile_20):
    """
    Fold80 is a way to calculate the evenness of coverage. It expresses the amount of additional sequencing needed to have
    80% of all targets covered at the currently observed mean. It is computed as the mean coverage divided by the 20th
    percentile.

    INPUT
    - coverage file: Output of Picard Raw WGS metrics file
    - Mean coverage: Mean coverage, parsed from the Raw WGS metrics file
    - Percentile_20: 20% of the total bases (genome_territory)
    """
    with open(coverage_file, 'r') as handle:
        lines = handle.readlines()
        histogram_start_index = lines.index("## HISTOGRAM\tjava.lang.Integer\n")

    # Read the entire coverage histogram.
    ## columns: coverage - coverage_count - unfiltered baseq count
    coverage_histo = lines[histogram_start_index+3:-1]
    coverage_histo = [l.rstrip().split("\t") for l in coverage_histo]

    # get the list of coverage values. As we do not change the order, the index is the coverage
    # e.g. [7, 100, 500, 100] -> 7 bases 0X, 100 bases 1X, 500 bases 2X, 100 bases 1X
    cov_values = [int(l[1]) for l in coverage_histo]

    # calculate cumsum for the histogram:
    # e.g. [7, 100, 500, 100] -> [7, 107, 607, 707]: 7 bases 0X, 107 bases max. 1X, 607 bases max. 2X, 707 bases max. 2X
    histo_cumsum = [sum(cov_values[:i+1]) for i, value in enumerate(cov_values)]

    # add the input percentile_20 to the cumsum list and sort
    tmp = histo_cumsum + [int(percentile_20)]
    tmp.sort()
    # Find the coverage of the 20th percentile
    # e.g. 20% 707 = 141 -> [7, 107, 141, 607, 707]
    find_index = tmp.index(int(percentile_20))
    val_before = histo_cumsum[find_index - 1]
    val_after = histo_cumsum[find_index]
    ## in our example, find_index will be 2, val_before: 107, val_after: 607

    cov_20 = (((val_after - percentile_20) / (val_after - val_before)) + find_index - 1)

    # ideally, as close as possible to 1
    return mean_coverage / cov_20