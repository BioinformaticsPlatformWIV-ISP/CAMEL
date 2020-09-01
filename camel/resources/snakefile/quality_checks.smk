from pathlib import Path

from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import quality_checks, trimming, contamination_check_kraken


rule quality_checks_kraken:
    """
    Checks the kraken output for contaminants.
    """
    input:
        INFORMS = (Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_INFORMS) if 'kraken' in config['analyses'] else []
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'kraken.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['kraken']
    run:
        import json
        if len(input) == 0:
            data_export = params.qc_check.to_dict()
        else:
            informs = SnakemakeUtils.load_object(input.INFORMS)
            contaminants = informs['contaminants_warn'] + informs['contaminants_fail']
            max_level = max([float(x[-1]) for x in contaminants] + [0.0])
            data_export = params.qc_check.to_dict(max_level)
        with open(output.JSON, 'w') as handle:
            json.dump(data_export, handle)


rule quality_checks_mapping_rate:
    """
    Checks the mapping rate against the reference genome / assembled contigs.
    """
    input:
        INFORMS = quality_checks.get_mapping_rate_informs(config)
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'mapping_ref.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY[f'map_rate_{config["quality_checks"].get("coverage_mode", "assembly")}']
    run:
        import json
        informs_bt2 = SnakemakeUtils.load_object(input.INFORMS)
        mapping_rate = float(informs_bt2['stats_map_rate'])
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(mapping_rate), handle)


rule quality_checks_coverage:
    """
    Checks the coverage against the reference genome / assembled contigs.
    """
    input:
        INFORMS = quality_checks.get_depth_informs(config)
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'cov_ref.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY[f'cov_{config["quality_checks"].get("coverage_mode", "assembly")}'],
        running_dir = Path(config['working_dir']) / 'quality_checks'
    run:
        import json
        samtools_depth_informs = SnakemakeUtils.load_object(input.INFORMS)
        median_depth = samtools_depth_informs['median_depth']
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(median_depth), handle)


rule quality_checks_typing_loci:
    """
    Checks the fraction of detected typing loci.
    """
    input:
        INFORMS = str(Path(config['working_dir']) / 'typing' / '{scheme}' / 'stats' / 'informs.io').format(
            scheme=config['quality_checks']['typing_scheme'])
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'cgmlst.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['cgmlst']
    run:
        import json
        cgmlst_stats = SnakemakeUtils.load_object(input.INFORMS)
        fraction_detected = cgmlst_stats['hits_found'] / cgmlst_stats['nb_of_loci']
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(100 * fraction_detected), handle)


rule quality_checks_parse_fastqc:
    """
    Tests additional quality metrics based on the FastQC data file output.
    """
    input:
        TXT = trimming.get_trimming_fastqc('post', config),
        TXT_RAW = trimming.get_trimming_fastqc('pre', config)
    output:
        INFORMS = Path(config['working_dir']) / 'quality_checks' / 'parse_fastqc' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'quality_checks' / 'parse_fastqc'
    run:
        from camel.app.tools.fastqc.fastqcdatafileparser import FastQCDataFileParser
        fastqc_checks = FastQCDataFileParser(camel)
        step = Step(rule, fastqc_checks, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(fastqc_checks, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(fastqc_checks, output)


rule quality_checks_fqc_gc_content:
    """
    Checks the additional QC checks based on the FastQC output.
    """
    input:
        INFORMS = rules.quality_checks_parse_fastqc.output.INFORMS
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'fastqc_gc_{ori}.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['fqc_gc'],
        gc_content_ref = config['quality_checks']['expected_gc_content'],
        index = lambda wildcards: 0 if (wildcards.ori == 'fwd') else 1
    run:
        import json
        informs = SnakemakeUtils.load_object(input.INFORMS)
        gc_diff = abs(params.gc_content_ref - informs['by_file']['gc_content'][params.index])
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(gc_diff), handle)


rule quality_checks_fqc_avg_read_quality:
    """
    Checks the additional QC checks based on the FastQC output.
    """
    input:
        INFORMS = rules.quality_checks_parse_fastqc.output.INFORMS
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'fastqc_avg_qual_{ori}.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['fqc_avg_qual'],
        index = lambda wildcards: 0 if (wildcards.ori == 'fwd') else 1
    run:
        import json
        informs = SnakemakeUtils.load_object(input.INFORMS)
        avg_quality = informs['by_file']['avg_read_qual'][params.index]
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(avg_quality), handle)


rule quality_checks_fqc_max_n_fraction:
    """
    Checks the maximal N-fraction based on the FastQC output.
    """
    input:
        INFORMS = rules.quality_checks_parse_fastqc.output.INFORMS
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'fastqc_n_fraction_{ori}.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['fqc_n_fraction'],
        index = lambda wildcards: 0 if (wildcards.ori == 'fwd') else 1
    run:
        import json
        informs = SnakemakeUtils.load_object(input.INFORMS)
        max_n_fraction = informs['by_file']['max_n_frac'][params.index]
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(max_n_fraction), handle)


rule quality_checks_fqc_seq_len:
    """
    Checks the sequence length distribution based on the FastQC output.
    """
    input:
        INFORMS = rules.quality_checks_parse_fastqc.output.INFORMS
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'fastqc_seq_len_{ori}.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['fqc_seq_len'],
        index = lambda wildcards: 0 if (wildcards.ori == 'fwd') else 1
    run:
        import json
        informs = SnakemakeUtils.load_object(input.INFORMS)
        perc = 100 * informs['by_file']['median_seq_len'][params.index] / informs['stats']['mode_read_length_raw']
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(perc), handle)


rule quality_checks_fqc_per_base:
    """
    Checks the per base sequence content based on the FastQC output.
    """
    input:
        INFORMS = rules.quality_checks_parse_fastqc.output.INFORMS
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'fastqc_per_base_{ori}.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['fqc_per_base'],
        index = lambda wildcards: 0 if (wildcards.ori == 'fwd') else 1
    run:
        import json
        informs = SnakemakeUtils.load_object(input.INFORMS)
        max_diff = informs['by_file']['max_per_base_diff'][params.index]
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(max_diff), handle)


rule quality_checks_fqc_qscore:
    """
    Checks the Q-score drop based on the FastQC output.
    """
    input:
        INFORMS = rules.quality_checks_parse_fastqc.output.INFORMS
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'fastqc_qscore_{ori}.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['fqc_qscore'],
        index = lambda wildcards: 0 if (wildcards.ori == 'fwd') else 1
    run:
        import json
        informs = SnakemakeUtils.load_object(input.INFORMS)
        qscore_drop = informs['by_file']['qscore_drop_pos'][params.index]
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(qscore_drop), handle)


rule quality_checks_report:
    """
    Creates the report for the quality checks workflow.
    """
    input:
        INFORMS = rules.quality_checks_parse_fastqc.output.INFORMS,
        JSON_kraken = rules.quality_checks_kraken.output.JSON if 'kraken' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_cgmlst = rules.quality_checks_typing_loci.output.JSON if 'typing' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_cov_ref = rules.quality_checks_coverage.output.JSON if 'coverage' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_map_rate_ref = rules.quality_checks_mapping_rate.output.JSON if 'coverage' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_avg_qual_fwd = rules.quality_checks_fqc_avg_read_quality.output.JSON.format(ori='fwd') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_avg_qual_rev = rules.quality_checks_fqc_avg_read_quality.output.JSON.format(ori='rev') if (config.get('read_type', 'illumina') == 'illumina' and 'fastqc' not in config['quality_checks'].get('disabled_checks', [])) else [],
        JSON_fqc_gc_fwd = rules.quality_checks_fqc_gc_content.output.JSON.format(ori='fwd') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_gc_rev = rules.quality_checks_fqc_gc_content.output.JSON.format(ori='rev') if (config.get('read_type', 'illumina') == 'illumina' and 'fastqc' not in config['quality_checks'].get('disabled_checks', [])) else [],
        JSON_fqc_n_frac_fwd = rules.quality_checks_fqc_max_n_fraction.output.JSON.format(ori='fwd') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_n_frac_rev = rules.quality_checks_fqc_max_n_fraction.output.JSON.format(ori='rev') if (config.get('read_type', 'illumina') == 'illumina' and 'fastqc' not in config['quality_checks'].get('disabled_checks', [])) else [],
        JSON_fqc_per_base_fwd = rules.quality_checks_fqc_per_base.output.JSON.format(ori='fwd') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_per_base_rev = rules.quality_checks_fqc_per_base.output.JSON.format(ori='rev') if (config.get('read_type', 'illumina') == 'illumina' and 'fastqc' not in config['quality_checks'].get('disabled_checks', [])) else [],
        JSON_fqc_qscore_fwd = rules.quality_checks_fqc_qscore.output.JSON.format(ori='fwd') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_qscore_rev = rules.quality_checks_fqc_qscore.output.JSON.format(ori='rev') if (config.get('read_type', 'illumina') == 'illumina' and 'fastqc' not in config['quality_checks'].get('disabled_checks', [])) else [],
        JSON_fqc_seq_len_fwd = rules.quality_checks_fqc_seq_len.output.JSON.format(ori='fwd') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_seq_len_rev = rules.quality_checks_fqc_seq_len.output.JSON.format(ori='rev') if (config.get('read_type', 'illumina') == 'illumina' and 'fastqc' not in config['quality_checks'].get('disabled_checks', [])) else []
    output:
        VAL_HTML = Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT,
        JSON = Path(config['working_dir']) / 'report' / 'informs.json'
    params:
        gc_content_ref = config['quality_checks']['expected_gc_content']
    run:
        # Combine the informs
        informs = []
        for json_file in input[1:]:
            with open(json_file) as handle:
                data = json.load(handle)
                data['ori'] = 'fwd' if 'fwd' in json_file else ('rev' if 'rev' in json_file else None)
                informs.append(data)
        with open(output.JSON, 'w') as handle:
            json.dump(informs, handle)

        # Create the report
        from camel.app.tools.pipelines.quality_checks.htmlreporterqualitychecks import HtmlReporterQualityChecks
        reporter = HtmlReporterQualityChecks(camel)
        reporter.update_parameters(gc_content_ref=params.gc_content_ref)
        reporter.add_input_informs({'qc_checks': informs, 'fastqc_parser': SnakemakeUtils.load_object(input.INFORMS)})
        reporter.run(str(Path(output.VAL_HTML).parent))
        SnakemakeUtils.dump_tool_output(reporter, 'VAL_HTML', output.VAL_HTML)


rule quality_checks_export_summary_info:
    """
    Exports the summary information for the quality checks workflow.
    """
    input:
        JSON = rules.quality_checks_report.output.JSON
    output:
        TSV = Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY
    run:
        with open(input.JSON) as handle:
            informs = json.load(handle)

        with open(output.TSV, 'w') as handle:
            for qc_check_data in informs:
                if qc_check_data['ori'] is None:
                    basename = f"qc_{qc_check_data['key']}"
                else:
                    basename = f"qc_{qc_check_data['key']}_{qc_check_data['ori']}"
                handle.write('\t'.join([f"{basename}_status", qc_check_data['status']]))
                handle.write('\n')
                handle.write('\t'.join([
                    f"{basename}_value", qc_check_data['fmt_string_value'].format(qc_check_data['value']) if
                        qc_check_data['value'] is not None else 'NA']))
                handle.write('\n')
