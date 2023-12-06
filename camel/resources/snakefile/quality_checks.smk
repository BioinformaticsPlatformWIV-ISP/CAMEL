from pathlib import Path
from typing import List, Union

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import quality_checks, contamination_check_kraken, trimming_ont, \
    quast, confindr, trimming_illumina


def _get_kraken2_informs(config, tech) -> Union[Path, List]:
    """
    Returns the informs used for the Kraken 2 informs.
    :param config: Snakemake configuration
    :param tech: Technology wildcards
    :return: Path to informs, or empty list if Kraken2 is disabled
    """
    # Kraken2 is disabled -> return empty list
    if 'kraken2' not in config['analyses']:
        return []

    # Illumina
    if tech == 'illumina':
        return Path(config['working_dir'], str(contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_INFORMS).format(
            read_key='fastq_pe'))
    # ONT
    elif tech == 'ont':
        return Path(config['working_dir'], str(contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_INFORMS).format(
            read_key='fastq_se'))
    else:
        raise ValueError(f"Invalid 'tech' wildcard: {tech}")

rule quality_checks_kraken:
    """
    Checks the kraken output for contaminants.
    """
    input:
        INFORMS = lambda wildcards: _get_kraken2_informs(config, wildcards.tech)
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'kraken_{tech}.json'
    params:
        qc_check = lambda wildcards: quality_checks.QC_CHECKS_BY_KEY[f'kraken_{wildcards.tech}']
    run:
        import json

        # noinspection PyTypeChecker
        if len(input) == 0:
            # noinspection PyUnresolvedReferences
            data_export = params.qc_check.to_dict()
        else:
            # noinspection PyTypeChecker
            informs = SnakemakeUtils.load_object(Path(input.INFORMS))
            contaminants = informs['contaminants_warn'] + informs['contaminants_fail']
            max_level = max([float(x[-1]) for x in contaminants] + [0.0])
            # noinspection PyUnresolvedReferences
            data_export = params.qc_check.to_dict(max_level)
        with open(output.JSON, 'w') as handle:
            json.dump(data_export, handle, indent=2)

rule quality_checks_typing_loci:
    """
    Checks the percentage of typing loci detected.
    """
    input:
        INFORMS = str(Path(config['working_dir']) / 'typing' / '{scheme}' / 'stats' / 'informs.io').format(
            scheme=config['quality_checks']['typing_scheme']) if config.get('quality_checks', {}).get('typing_scheme') is not None else []
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'typing_loci.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['typing_loci']
    run:
        import json

        if len(input.INFORMS) == 0:
            data_export = params.qc_check.to_dict()
        else:
            typing_stats = SnakemakeUtils.load_object(Path(input.INFORMS))
            fraction_detected = typing_stats['hits_found'] / typing_stats['nb_of_loci']
            data_export = params.qc_check.to_dict(100 * fraction_detected)
        with open(output.JSON, 'w') as handle:
            json.dump(data_export, handle, indent=2)

rule quality_checks_mapping_rate_se:
    """
    Checks the mapping rate against the reference genome / assembled contigs for the SE reads (if available).
    """
    input:
        INFORMS = quality_checks.get_mapping_rate_informs(config, 'fastq_se')
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'map_rate_{mode}_ont.json'
    params:
        qc_check = lambda wildcards: quality_checks.QC_CHECKS_BY_KEY[f'map_rate_{wildcards.mode}_ont'],
        key = 'mapping_perc'
    run:
        import json
        informs_mapping = SnakemakeUtils.load_object(Path(input.INFORMS))
        mapping_rate = float(informs_mapping[params.key])
        with open(output.JSON, 'w') as handle:
            # noinspection PyUnresolvedReferences
            json.dump(params.qc_check.to_dict(mapping_rate), handle, indent=2)

rule quality_checks_mapping_rate_pe:
    """
    Checks the mapping rate against the reference genome / assembled contigs for the PE reads (if available).
    """
    input:
        INFORMS = quality_checks.get_mapping_rate_informs(config, 'fastq_pe')
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'map_rate_{mode}_illumina.json'
    params:
        qc_check = lambda wildcards: quality_checks.QC_CHECKS_BY_KEY[f'map_rate_{wildcards.mode}_illumina'],
        key = 'stats_map_rate'
    run:
        import json
        informs_mapping = SnakemakeUtils.load_object(Path(input.INFORMS))
        mapping_rate = float(informs_mapping[params.key])
        with open(output.JSON, 'w') as handle:
            # noinspection PyUnresolvedReferences
            json.dump(params.qc_check.to_dict(mapping_rate), handle, indent=2)

rule quality_checks_depth_se:
    """
    Checks the coverage against the assembled contigs.
    """
    input:
        INFORMS = lambda wildcards: quality_checks.get_depth_informs(config, 'fastq_se')
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'cov_{mode}_ont.json'
    params:
        qc_check = lambda wildcards: quality_checks.QC_CHECKS_BY_KEY[f'cov_{wildcards.mode}_ont'],
        running_dir = Path(config['working_dir']) / 'quality_checks'
    run:
        import json

        # noinspection PyTypeChecker
        samtools_depth_informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        median_depth = samtools_depth_informs['median_depth']
        with open(output.JSON, 'w') as handle:
            # noinspection PyUnresolvedReferences
            json.dump(params.qc_check.to_dict(median_depth), handle, indent=2)

rule quality_checks_depth_pe:
    """
    Checks the coverage against the assembled contigs.
    """
    input:
        INFORMS = lambda wildcards: quality_checks.get_depth_informs(config, 'fastq_pe')
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'cov_{mode}_illumina.json'
    params:
        qc_check = lambda wildcards: quality_checks.QC_CHECKS_BY_KEY[f'cov_{wildcards.mode}_illumina'],
        running_dir = Path(config['working_dir']) / 'quality_checks'
    run:
        import json

        # noinspection PyTypeChecker
        samtools_depth_informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        import pprint
        pprint.pprint(samtools_depth_informs)
        median_depth = samtools_depth_informs['median_depth']
        with open(output.JSON, 'w') as handle:
            # noinspection PyUnresolvedReferences
            json.dump(params.qc_check.to_dict(median_depth), handle, indent=2)

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
        import json
        import pandas as pd

        data_quast = pd.read_table(input.TSV, names=['key', 'value'])
        if data_quast[data_quast['key'] == 'assembly_total_length_ref'].iloc[0]['value'] != '-':
            # Calc the percentage deviation from the reference genome length
            total_length = int(data_quast[data_quast['key'] == 'assembly_total_length'].iloc[0]['value'])
            total_length_ref = int(data_quast[data_quast['key'] == 'assembly_total_length_ref'].iloc[0]['value'])
            perc_deviation = abs(100 * ((total_length / total_length_ref) - 1))
        else:
            # Skip test if reference genome is not available
            perc_deviation = None
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(perc_deviation), handle, indent=2)

rule quality_checks_parse_fastqc:
    """
    Tests additional quality metrics based on the FastQC data file output.
    """
    input:
        TXT = Path(config['working_dir'], trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_FASTQC_TXT_POST),
        TXT_RAW = Path(config['working_dir'], trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_FASTQC_TXT_PRE)
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
        JSON = Path(config['working_dir']) / 'quality_checks' / 'fqc_gc_{ori}.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['fqc_gc_{ori}'],
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
        JSON = Path(config['working_dir']) / 'quality_checks' / 'fqc_avg_qual_{ori}.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['fqc_avg_qual_{ori}'],
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
        JSON = Path(config['working_dir']) / 'quality_checks' / 'fqc_n_fraction_{ori}.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['fqc_n_fraction_{ori}'],
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
        JSON = Path(config['working_dir']) / 'quality_checks' / 'fqc_seq_len_{ori}.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['fqc_seq_len_{ori}'],
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
        JSON = Path(config['working_dir']) / 'quality_checks' / 'fqc_per_base_{ori}.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['fqc_per_base_{ori}'],
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
        JSON = Path(config['working_dir']) / 'quality_checks' / 'fqc_qscore_{ori}.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['fqc_qscore_{ori}'],
        index = lambda wildcards: 0 if (wildcards.ori == 'fwd') else 1
    run:
        import json
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        qscore_drop = informs['by_file']['qscore_drop_pos'][params.index]
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(qscore_drop), handle, indent=2)

rule quality_checks_confindr:
    """
    Extracts the quality checks from the ConFindr output.
    """
    input:
        INFORMS = Path(config['working_dir']) / confindr.OUTPUT_CONFINDR_INFORMS if 'confindr' in config['analyses'] else []
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'confindr.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['confindr']
    run:
        import json
        if len(input.INFORMS) == 0:
            data_export = params.qc_check.to_dict()
        else:
            informs = SnakemakeUtils.load_object(Path(input.INFORMS))
            nb_contam_snps = informs['NumContamSNVs']
            data_export = params.qc_check.to_dict(nb_contam_snps)
        with open(output.JSON, 'w') as handle:
            json.dump(data_export, handle, indent=2)

rule quality_checks_busco:
    """
    Extracts the BUSCO metric from the QUAST output.
    """
    input:
        INFORMS = Path(config['working_dir']) / quast.OUTPUT_BUSCO_INFORMS
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'busco.json'
    params:
        qc_check = quality_checks.QC_CHECKS_BY_KEY['busco']
    run:
        import json
        data_in = SnakemakeUtils.load_object(Path(input.INFORMS))
        value = data_in['results']['results']['Complete']
        with open(output.JSON, 'w') as handle:
            json.dump(params.qc_check.to_dict(value), handle, indent=2)

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
        JSON = Path(config['working_dir']) / 'quality_checks' / 'seqkit_gc.json'
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

rule quality_checks_combine_all:
    """
    Collects the quality checks for read type 'illumina'.
    """
    input:
        JSON = [Path(config['working_dir'], path_json) for path_json in quality_checks.get_qc_checks(
            config['input_type'], config.get('quality_checks', {}).get('skipped', []))]
    output:
        JSON = Path(config['working_dir']) / 'quality_checks' / 'all.json'
    run:
        import json

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

rule quality_checks_report:
    """
    Creates the report for the quality checks workflow.
    """
    input:
        JSON = rules.quality_checks_combine_all.output.JSON
    output:
        VAL_HTML = Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT
    run:
        import json

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
        JSON =  rules.quality_checks_combine_all.output.JSON
    output:
        TSV = Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY
    run:
        import json

        with open(input.JSON) as handle:
            informs = json.load(handle)

        with open(output.TSV, 'w') as handle:
            for qc_check_data in informs:
                if qc_check_data.get('ori') is None:
                    basename = f"qc_{qc_check_data['key']}"
                else:
                    basename = f"qc_{qc_check_data['key'].format(ori=qc_check_data['ori'])}"
                handle.write('\t'.join([f"{basename}_status", qc_check_data['status']]))
                handle.write('\n')
                handle.write('\t'.join([
                    f"{basename}_value", qc_check_data['fmt_string_value'].format(qc_check_data['value']) if
                        qc_check_data['value'] is not None else 'NA']))
                handle.write('\n')
