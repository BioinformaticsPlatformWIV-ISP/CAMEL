from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import quality_checks, trimming, contamination_check_kraken, trimming_ont, \
    assembly_spades, quast


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

        # noinspection PyTypeChecker
        if len(input) == 0:
            data_export = params.qc_check.to_dict()
        else:
            informs = SnakemakeUtils.load_object(Path(input.INFORMS))
            contaminants = informs['contaminants_warn'] + informs['contaminants_fail']
            max_level = max([float(x[-1]) for x in contaminants] + [0.0])
            data_export = params.qc_check.to_dict(max_level)
        with open(output.JSON, 'w') as handle:
            json.dump(data_export, handle, indent=2)

rule quality_checks_mapping_rate:
    """
    Checks the mapping rate against the reference genome / assembled contigs.
    """
    input:
        INFORMS = quality_checks.get_mapping_rate_informs(config)
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'mapping_ref.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY[f'map_rate_{config["quality_checks"].get("coverage_mode", "assembly")}'],
        key = 'stats_map_rate' if config['read_type'] == 'illumina' else 'mapping_perc'
    run:
        import json
        informs_mapping = SnakemakeUtils.load_object(Path(input.INFORMS))
        mapping_rate = float(informs_mapping[params.key])
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(mapping_rate), handle, indent=2)

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
        samtools_depth_informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        median_depth = samtools_depth_informs['median_depth']
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(median_depth), handle, indent=2)

rule quality_checks_typing_loci:
    """
    Checks the fraction of detected typing loci.
    """
    input:
        INFORMS = str(Path(config['working_dir']) / 'typing' / '{scheme}' / 'stats' / 'informs.io').format(
            scheme=config['quality_checks']['typing_scheme']) if 'typing' not in config['quality_checks'].get('disabled_checks', []) else []
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'cgmlst.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['cgmlst']
    run:
        import json
        cgmlst_stats = SnakemakeUtils.load_object(Path(input.INFORMS))
        fraction_detected = cgmlst_stats['hits_found'] / cgmlst_stats['nb_of_loci']
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(100 * fraction_detected), handle, indent=2)

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
        fastqc_checks = FastQCDataFileParser(Camel.get_instance())
        step = Step(str(rule), fastqc_checks, Camel.get_instance(), params.running_dir)
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
        import dataclasses
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        gc_diff = abs(params.gc_content_ref - informs['by_file']['gc_content'][params.index])

        # Fill in expected GC content in explanation
        qc_check = dataclasses.replace(params.qc_check,
            explanation=params.qc_check.explanation.format(params.gc_content_ref))

        with open(output.JSON, 'w') as handle:
            json.dump(qc_check.to_dict(gc_diff), handle, indent=2)

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
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        avg_quality = informs['by_file']['avg_read_qual'][params.index]
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(avg_quality), handle, indent=2)

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
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        max_n_fraction = informs['by_file']['max_n_frac'][params.index]
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(max_n_fraction), handle, indent=2)

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
        import dataclasses

        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        mode_length = informs['stats']['mode_read_length_raw']
        perc = 100 * informs['by_file']['median_seq_len'][params.index] / mode_length

        # Fill in value in parameter explanation
        qc_check = dataclasses.replace(params.qc_check, explanation=params.qc_check.explanation.format(mode_length))

        # Export as JSON
        with open(output.JSON, 'w') as handle:
            json.dump(qc_check.to_dict(perc), handle, indent=2)

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
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        max_diff = informs['by_file']['max_per_base_diff'][params.index]
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(max_diff), handle, indent=2)

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
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        qscore_drop = informs['by_file']['qscore_drop_pos'][params.index]
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(qscore_drop), handle, indent=2)

rule quality_checks_seqkit_stats:
    """
    Runs seqkit stats to extract statistics for Nanopore data.
    """
    input:
        FASTQ = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_READS
    output:
        INFORMS = Path(config['working_dir']) / 'quality_checks' / 'seqkit' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'quality_checks' / 'seqkit'
    run:
        from camel.app.tools.seqkit.seqkitstats import SeqkitStats
        seqkit_stats = SeqkitStats(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(seqkit_stats, input)
        step = Step(str(rule), seqkit_stats, Camel.get_instance(), params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqkit_stats, output)

rule quality_checks_seqkit_gc:
    """
    Checks the %GC-content on the seqkit output.
    """
    input:
        INFORMS = rules.quality_checks_seqkit_stats.output.INFORMS
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'nanoplot_gc.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['seqkit_gc'],
        gc_content_ref = config['quality_checks']['expected_gc_content']
    run:
        import json
        import dataclasses
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        gc_diff = abs(params.gc_content_ref - informs['GC(%)'])

        # Fill in value in parameter explanation
        qc_check = dataclasses.replace(params.qc_check, explanation=params.qc_check.explanation.format(
            params.gc_content_ref))

        with open(output.JSON, 'w') as handle:
            json.dump(qc_check.to_dict(gc_diff), handle, indent=2)

rule quality_checks_nanoplot_qual:
    """
    Checks the read quality based on the NanoPlot output.
    """
    input:
        INFORMS = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_INFORMS_POST
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'nanoplot_qual.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['nanoplot_qual']
    run:
        import json
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(float(informs['median_qual'])), handle, indent=2)

rule quality_checks_nanoplot_len:
    """
    Checks the read length based on the NanoPlot output.
    """
    input:
        INFORMS = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_NANOPLOT_INFORMS_POST
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'nanoplot_len.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['nanoplot_len']
    run:
        import json
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(int(float(informs['median_read_length']))), handle, indent=2)

rule quality_checks_assembly_total_len:
    """
    Checks if the total assembly length is within the expected range.
    """
    input:
        TSV = Path(config['working_dir']) / quast.OUTPUT_QUAST_SUMMARY
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'assembly_total_len.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['assembly_total_len']
    run:
        import pandas as pd

        # Calc the percentage deviation from the reference genome length
        data_quast = pd.read_table(input.TSV, names=['key', 'value'])
        total_length = data_quast[data_quast['key'] == 'assembly_total_length'].iloc[0]['value']
        total_length_ref = data_quast[data_quast['key'] == 'assembly_total_length_ref'].iloc[0]['value']
        perc_deviation = abs(100 * ((total_length / total_length_ref) - 1))

        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(perc_deviation), handle, indent=2)

rule quality_checks_combine_illumina:
    """
    Collects the quality checks for read type 'illumina'.
    """
    input:
        JSON_kraken = rules.quality_checks_kraken.output.JSON if 'kraken' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_cgmlst = rules.quality_checks_typing_loci.output.JSON if 'typing' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_cov_ref = rules.quality_checks_coverage.output.JSON if 'coverage' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_map_rate_ref = rules.quality_checks_mapping_rate.output.JSON if 'coverage' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_assembly_total_len = rules.quality_checks_assembly_total_len.output.JSON if 'assembly' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_avg_qual_fwd = str(rules.quality_checks_fqc_avg_read_quality.output.JSON).format(ori='fwd') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else[],
        JSON_fqc_avg_qual_rev = str(rules.quality_checks_fqc_avg_read_quality.output.JSON).format(ori='rev') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_gc_fwd = str(rules.quality_checks_fqc_gc_content.output.JSON).format(ori='fwd') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_gc_rev = str(rules.quality_checks_fqc_gc_content.output.JSON).format(ori='rev') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_n_frac_fwd = str(rules.quality_checks_fqc_max_n_fraction.output.JSON).format(ori='fwd') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_n_frac_rev = str(rules.quality_checks_fqc_max_n_fraction.output.JSON).format(ori='rev') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_per_base_fwd = str(rules.quality_checks_fqc_per_base.output.JSON).format(ori='fwd') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_per_base_rev = str(rules.quality_checks_fqc_per_base.output.JSON).format(ori='rev') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_qscore_fwd = str(rules.quality_checks_fqc_qscore.output.JSON).format(ori='fwd') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_qscore_rev = str(rules.quality_checks_fqc_qscore.output.JSON).format(ori='rev') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_seq_len_fwd = str(rules.quality_checks_fqc_seq_len.output.JSON).format(ori='fwd') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_fqc_seq_len_rev = str(rules.quality_checks_fqc_seq_len.output.JSON).format(ori='rev') if 'fastqc' not in config['quality_checks'].get('disabled_checks', []) else []
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'illumina.json'
    run:
        # Combine QC check informs
        informs_out = []
        for path_json in [Path(x) for x in input]:
            with path_json.open() as handle:
                data = json.load(handle)
            data['ori'] = 'fwd' if 'fwd' in path_json.name else ('rev' if 'rev' in path_json.name else None)
            informs_out.append(data)

        # Save output file
        with open(output.JSON, 'w') as handle:
            json.dump(informs_out, handle, indent=2)

rule quality_checks_combine_nanopore:
    """
    Collects the quality checks for read type 'nanopore'.
    """
    input:
        JSON_kraken = rules.quality_checks_kraken.output.JSON if 'kraken' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_cgmlst = rules.quality_checks_typing_loci.output.JSON if 'typing' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_cov_ref = rules.quality_checks_coverage.output.JSON if 'coverage' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_map_rate_ref = rules.quality_checks_mapping_rate.output.JSON if 'coverage' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_nanoplot_len = rules.quality_checks_nanoplot_len.output.JSON if 'length' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_nanoplot_qual = rules.quality_checks_nanoplot_qual.output.JSON if 'quality' not in config['quality_checks'].get('disabled_checks', []) else [],
        JSON_nanoplot_gc = rules.quality_checks_seqkit_gc.output.JSON if 'gc' not in config['quality_checks'].get('disabled_checks', []) else []
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'nanopore.json'
    run:
        # Combine QC check informs
        informs_out = []
        for path_json in [Path(x) for x in input]:
            with path_json.open() as handle:
                data = json.load(handle)
            informs_out.append(data)

        # Save output file
        with open(output.JSON, 'w') as handle:
            json.dump(informs_out, handle, indent=2)

rule quality_checks_report:
    """
    Creates the report for the quality checks workflow.
    """
    input:
        JSON = Path(config['working_dir']) / 'quality_checks' / f"{config['read_type']}.json"
    output:
        VAL_HTML = Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT
    run:
        from camel.app.tools.pipelines.quality_checks.htmlreporterqualitychecks import HtmlReporterQualityChecks
        reporter = HtmlReporterQualityChecks(Camel.get_instance())
        with open(input.JSON) as handle:
            informs = json.load(handle)
        reporter.add_input_informs({'qc_checks': informs})
        reporter.run(Path(output.VAL_HTML).parent)
        SnakemakeUtils.dump_tool_output(reporter, 'VAL_HTML', Path(output.VAL_HTML))

rule quality_checks_export_summary_info:
    """
    Exports the summary information for the quality checks workflow.
    """
    input:
        JSON = Path(config['working_dir']) / 'quality_checks' / f"{config['read_type']}.json"
    output:
        TSV = Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY
    run:
        with open(input.JSON) as handle:
            informs = json.load(handle)

        with open(output.TSV, 'w') as handle:
            for qc_check_data in informs:
                if qc_check_data.get('ori') is None:
                    basename = f"qc_{qc_check_data['key']}"
                else:
                    basename = f"qc_{qc_check_data['key']}_{qc_check_data['ori']}"
                handle.write('\t'.join([f"{basename}_status", qc_check_data['status']]))
                handle.write('\n')
                handle.write('\t'.join([
                    f"{basename}_value", qc_check_data['fmt_string_value'].format(qc_check_data['value']) if
                        qc_check_data['value'] is not None else 'NA']))
                handle.write('\n')
