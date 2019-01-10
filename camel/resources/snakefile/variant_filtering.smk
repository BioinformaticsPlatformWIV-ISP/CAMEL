import logging

import os
import shutil

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.variant_calling import OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ
from camel.resources.snakefile.variant_filtering import OUTPUT_VARIANT_FILTERING_STATS, OUTPUT_VARIANT_FILTERING_VCF, \
    OUTPUT_VARIANT_FILTERING_SUMMARY, get_filtering_param

camel = Camel.get_instance()


def filter_is_disabled(filter_key: str) -> bool:
    """
    Returns True if the given filter is disabled, False otherwise.
    :param filter_key: Filter key
    :return: True if disabled, false otherwise
    """
    if filter_key not in config['variant_filtering']:
        return False
    return config['variant_filtering'][filter_key].get('disabled', False)


rule Variant_filtering_depth:
    """
    Filters variants based on sequencing depth.
    """
    input:
        VCF_GZ=os.path.join(config['working_dir'], OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ)
    output:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'depth', 'vcf_gz.io'),
        INFORMS=os.path.join(config['working_dir'], 'variant_filtering', 'depth', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_filtering', 'depth'),
        min_total_depth=get_filtering_param(config, 'depth', 'min_total_depth'),
        min_fwd_depth=get_filtering_param(config, 'depth', 'min_fwd_depth'),
        min_rev_depth=get_filtering_param(config, 'depth', 'min_rev_depth')
    run:
        from camel.app.tools.variantfiltering.depthfilter import DepthFilter
        depth_filter = DepthFilter(camel)
        step = Step(rule, depth_filter, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(depth_filter, input)
        if params.min_total_depth is not None:
            depth_filter.update_parameters(min_depth=params.min_total_depth)
        if params.min_fwd_depth is not None:
            depth_filter.update_parameters(min_forward_depth=params.min_fwd_depth)
        if params.min_rev_depth is not None:
            depth_filter.update_parameters(min_reverse_depth=params.min_rev_depth)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(depth_filter, output)

rule Variant_filtering_snp_quality:
    """
    Filters variants based on SNP quality.
    """
    input:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'depth', 'vcf_gz.io'),
    output:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'snp_quality', 'vcf_gz.io'),
        INFORMS=os.path.join(config['working_dir'], 'variant_filtering', 'snp_quality', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_filtering', 'snp_quality'),
        min_snp_quality=get_filtering_param(config, 'snp_quality', 'min_snp_quality')
    run:
        from camel.app.tools.variantfiltering.snpqualityfilter import SnpQualityFilter
        snp_qual_filter = SnpQualityFilter(camel)
        SnakemakeUtils.add_pickle_inputs(snp_qual_filter, input)
        step = Step(rule, snp_qual_filter, camel, params.running_dir, config)
        if params.min_snp_quality is not None:
            snp_qual_filter.update_parameters(min_snp_quality=params.min_snp_quality)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(snp_qual_filter, output)

rule Variant_filtering_mapping_quality:
    """
    Filters variants based on mapping quality.
    """
    input:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'snp_quality', 'vcf_gz.io')
    output:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'mapping_quality', 'vcf_gz.io'),
        INFORMS=os.path.join(config['working_dir'], 'variant_filtering', 'mapping_quality', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_filtering', 'mapping_quality'),
        min_mapping_quality=get_filtering_param(config, 'mapping_quality', 'min_mapping_quality')
    run:
        from camel.app.tools.variantfiltering.mappingqualityfilter import MappingQualityFilter
        snp_qual_filter = MappingQualityFilter(camel)
        SnakemakeUtils.add_pickle_inputs(snp_qual_filter, input)
        step = Step(rule, snp_qual_filter, camel, params.running_dir, config)
        if params.min_mapping_quality is not None:
            snp_qual_filter.update_parameters(min_mapping_quality=params.min_mapping_quality)
        snp_qual_filter.update_parameters(output_filename="filtered_mapping_qual.vcf.gz")
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(snp_qual_filter, output)

rule Variant_filtering_distance_index:
    """
    Indexes the VCF file before applying the distance filter.
    """
    input:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'mapping_quality', 'vcf_gz.io')
    output:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'distance', 'vcf_gz-indexed.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_filtering', 'distance')
    run:
        from camel.app.tools.bcftools.bcftoolsindex import BcftoolsIndex
        bcftools_index  = BcftoolsIndex(camel)
        step = Step(rule, bcftools_index, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(bcftools_index, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_index, output)

rule Variant_filtering_distance:
    """
    Filters the variants based on distance.
    """
    input:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'distance', 'vcf_gz-indexed.io')
    output:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'distance', 'vcf_gz.io'),
        INFORMS=os.path.join(config['working_dir'], 'variant_filtering', 'distance', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_filtering', 'distance'),
        min_distance=get_filtering_param(config, 'distance', 'min_distance'),
        keep_best=get_filtering_param(config, 'distance', 'keep_best')
    run:
        from camel.app.tools.variantfiltering.distancefilter import DistanceFilter
        distance_filter = DistanceFilter(camel)
        step = Step(rule, distance_filter, camel, params.running_dir, config)
        if params.min_distance is not None:
            distance_filter.update_parameters(min_distance=params.min_distance)
        if params.keep_best is not None:
            distance_filter.update_parameters(keep_best=params.keep_best)
        SnakemakeUtils.add_pickle_inputs(distance_filter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(distance_filter, output)

rule Variant_filtering_zscore_index:
    """
    Indexes the VCF file before applying the Z-score filter.
    """
    input:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'distance', 'vcf_gz.io')
    output:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'zscore', 'vcf_gz-indexed.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_filtering', 'zscore')
    run:
        from camel.app.tools.bcftools.bcftoolsindex import BcftoolsIndex
        bcftools_index  = BcftoolsIndex(camel)
        step = Step(rule, bcftools_index, camel, params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(bcftools_index, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_index, output)

rule Variant_filtering_zscore:
    """
    Filters variants based on the Z-score metric.
    """
    input:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'zscore', 'vcf_gz-indexed.io'),
        BAM=os.path.join(config['working_dir'], 'variant_calling', 'alignment_sorting', 'bam-sorted.io')
    output:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'zscore', 'vcf_gz.io'),
        INFORMS=os.path.join(config['working_dir'], 'variant_filtering', 'zscore', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_filtering', 'zscore'),
        min_zscore=get_filtering_param(config, 'zscore', 'min_zscore'),
        y_multiplier=get_filtering_param(config, 'zscore', 'y_multiplier')
    run:
        from camel.app.tools.variantfiltering.zscorefilter import ZScoreFilter
        zscore_filter = ZScoreFilter(camel)
        step = Step(rule, zscore_filter, camel, params.running_dir, config)
        if params.min_zscore is not None:
            zscore_filter.update_parameters(min_zscore=params.min_zscore)
        if params.y_multiplier is not None:
            zscore_filter.update_parameters(y_multiplier=params.y_multiplier)

        if len(SnakemakeUtils.load_object(input.BAM)) == 0:
            logging.info("No BAM input found, skipping Z-score filter.")
            shutil.copyfile(input.VCF_GZ, output.VCF_GZ)
            SnakemakeUtils.dump_object({'variants_in': 'NA', 'variants_out': 'NA'}, output.INFORMS)
        else:
            SnakemakeUtils.add_pickle_inputs(zscore_filter, input)
            step.run_step()
            SnakemakeUtils.dump_tool_outputs(zscore_filter, output)

rule Variant_filtering_unzip_vcf_zscore:
    """
    Unzips the filtered VCF file.
    """
    input:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'zscore', 'vcf_gz.io')
    output:
        VCF=os.path.join(config['working_dir'], OUTPUT_VARIANT_FILTERING_VCF)
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_filtering', 'unzip')
    run:
        from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
        bcftools_view = BcftoolsView(camel)
        SnakemakeUtils.add_pickle_inputs(bcftools_view, input)
        step = Step(rule, bcftools_view, camel, params.running_dir, config)
        bcftools_view.update_parameters(output_format='VCF', compress_output=False, output_filename='variants_filtered.vcf')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_view, output)

rule Variant_filtering_region:
    """
    Filters out regions that are problematic for variant calling based on an input BED file.
    """
    input:
        VCF_GZ=os.path.join(config['working_dir'], 'variant_filtering', 'zscore', 'vcf_gz.io'),
    output:
        VCF=os.path.join(config['working_dir'], 'variant_filtering', 'regions', 'vcf.io'),
        INFORMS=os.path.join(config['working_dir'], 'variant_filtering', 'regions', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'variant_filtering', 'regions'),
        bed_file=get_filtering_param(config, 'region', 'bed_file')
    run:
        from camel.app.components.vcf.vcfutils import VCFUtils
        from camel.app.tools.bcftools.bcftoolsfilter import BcftoolsFilter
        bcftools_filter = BcftoolsFilter(camel)
        vcf_file = SnakemakeUtils.load_object(input.VCF_GZ)[0]
        bcftools_filter.add_input_files({'VCF_GZ': [vcf_file]})
        if params.bed_file is not None:
            bcftools_filter.add_input_files({'BED': [ToolIOFile(params.bed_file)]})
        bcftools_filter.update_parameters(output_type='v', invert_targets=True,
                                          output_filename='variants_regions_removed.vcf')
        step = Step(rule, bcftools_filter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(bcftools_filter, 'VCF', output.VCF)
        informs = {
            'variants_out': VCFUtils.count_variants(bcftools_filter.tool_outputs['VCF'][0].path),
            'variants_in': VCFUtils.count_variants(vcf_file.path)
        }
        SnakemakeUtils.dump_object(informs, output.INFORMS)

rule Variant_filtering_collect_stats:
    """
    Collects the stats for all of the filtering steps.
    """
    input:
        INFORMS_depth=os.path.join(config['working_dir'], 'variant_filtering', 'depth', 'informs.io'),
        INFORMS_snp_qual=os.path.join(config['working_dir'], 'variant_filtering', 'snp_quality', 'informs.io'),
        INFORMS_mapping_qual=os.path.join(config['working_dir'], 'variant_filtering', 'mapping_quality', 'informs.io'),
        INFORMS_distance=os.path.join(config['working_dir'], 'variant_filtering', 'distance', 'informs.io'),
        INFORMS_zscore=os.path.join(config['working_dir'], 'variant_filtering', 'zscore', 'informs.io'),
        INFORMS_region=os.path.join(config['working_dir'], 'variant_filtering', 'regions', 'informs.io')
    output:
        JSON=os.path.join(config['working_dir'], OUTPUT_VARIANT_FILTERING_STATS)
    params:
        working_dir=os.path.join(config['working_dir'], 'variant_filtering', 'stats')
    run:
        import json
        filtering_data = {}
        for input_key in input.keys():
            filter_name = input_key.replace('INFORMS_', '')
            filtering_data[filter_name] = SnakemakeUtils.load_object(input[input_key])
        output_path = os.path.join(params.working_dir, 'stats.json')
        with open(output_path, 'w') as handle:
            json.dump(filtering_data, handle)
        SnakemakeUtils.dump_object([ToolIOFile(output_path)], output.JSON)

rule Variant_filtering_dump_summary_info:
    """
    Dumps the summary information for the variant filtering workflow.
    """
    input:
        JSON=os.path.join(config['working_dir'], OUTPUT_VARIANT_FILTERING_STATS)
    output:
        os.path.join(config['working_dir'], OUTPUT_VARIANT_FILTERING_SUMMARY)
    run:
        import json
        with open(SnakemakeUtils.load_object(input.JSON)[0].path) as handle:
            filtering_info = json.load(handle)

        with open(output[0], 'w') as handle:
            for key, data in sorted(filtering_info.items()):
                handle.write('\t'.join([f'filt-{key}-in', str(data['variants_in'])]))
                handle.write('\n')
                handle.write('\t'.join([f'filt-{key}-out', str(data['variants_out'])]))
                handle.write('\n')
