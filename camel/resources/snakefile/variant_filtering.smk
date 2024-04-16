from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import variant_calling, variant_filtering


rule variant_filtering_depth:
    """
    Filters variants based on sequencing depth.
    """
    input:
        VCF_GZ = Path(config['working_dir']) / variant_calling.get_vcf_gz(config)
    output:
        VCF_GZ = Path(config['working_dir']) / 'variant_filtering' / '01-depth' / 'vcf_gz.io',
        INFORMS = Path(config['working_dir']) / 'variant_filtering' / '01-depth'/ 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'variant_filtering' / '01-depth',
        min_total_depth = variant_filtering.get_filtering_param(config, 'depth', 'min_total_depth'),
        min_fwd_depth = variant_filtering.get_filtering_param(config, 'depth', 'min_fwd_depth'),
        min_rev_depth = variant_filtering.get_filtering_param(config, 'depth', 'min_rev_depth'),
        soft_filter = config['variant_filtering'].get('soft_filter', False)
    run:
        from camel.app.tools.variantfiltering.depthfilter import DepthFilter
        depth_filter = DepthFilter(Camel.get_instance())
        step = Step(str(rule), depth_filter, Camel.get_instance(), params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(depth_filter, input)
        if params.min_total_depth is not None:
            depth_filter.update_parameters(min_depth=params.min_total_depth)
        if params.min_fwd_depth is not None:
            depth_filter.update_parameters(min_forward_depth=params.min_fwd_depth)
        if params.min_rev_depth is not None:
            depth_filter.update_parameters(min_reverse_depth=params.min_rev_depth)
        depth_filter.update_parameters(soft_filter=params.soft_filter)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(depth_filter, output)

rule variant_filtering_snp_quality:
    """
    Filters variants based on SNP quality.
    """
    input:
        VCF_GZ = rules.variant_filtering_depth.output.VCF_GZ
    output:
        VCF_GZ = Path(config['working_dir']) / 'variant_filtering' / '02-snp_qual' / 'vcf_gz.io',
        INFORMS = Path(config['working_dir']) / 'variant_filtering' / '02-snp_qual' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'variant_filtering' / '02-snp_qual',
        min_snp_quality = variant_filtering.get_filtering_param(config, 'snp_quality', 'min_snp_quality'),
        soft_filter = config['variant_filtering'].get('soft_filter', False)
    run:
        from camel.app.tools.variantfiltering.snpqualityfilter import SnpQualityFilter
        snp_qual_filter = SnpQualityFilter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(snp_qual_filter, input)
        step = Step(str(rule), snp_qual_filter, Camel.get_instance(), params.running_dir, config)
        if params.min_snp_quality is not None:
            snp_qual_filter.update_parameters(min_snp_quality=params.min_snp_quality)
        snp_qual_filter.update_parameters(soft_filter=params.soft_filter)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(snp_qual_filter, output)

rule variant_filtering_mapping_quality:
    """
    Filters variants based on mapping quality.
    """
    input:
        VCF_GZ = rules.variant_filtering_snp_quality.output.VCF_GZ
    output:
        VCF_GZ = Path(config['working_dir']) / 'variant_filtering' / '03-mapping_qual' / 'vcf_gz.io',
        INFORMS = Path(config['working_dir']) / 'variant_filtering'/ '03-mapping_qual' /'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'variant_filtering' / '03-mapping_qual',
        min_mapping_quality = variant_filtering.get_filtering_param(config, 'mapping_quality', 'min_mapping_quality'),
        soft_filter = config['variant_filtering'].get('soft_filter', False)
    run:
        from camel.app.tools.variantfiltering.mappingqualityfilter import MappingQualityFilter
        mapping_quality_filter = MappingQualityFilter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(mapping_quality_filter, input)
        step = Step(str(rule), mapping_quality_filter, Camel.get_instance(), params.running_dir, config)
        if params.min_mapping_quality is not None:
            mapping_quality_filter.update_parameters(min_mapping_quality=params.min_mapping_quality)
        mapping_quality_filter.update_parameters(
            output_filename="filtered_mapping_qual.vcf.gz", soft_filter=params.soft_filter)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(mapping_quality_filter, output)

rule variant_filtering_distance_index:
    """
    Indexes the VCF file before applying the distance filter.
    """
    input:
        VCF_GZ = rules.variant_filtering_mapping_quality.output.VCF_GZ
    output:
        VCF_GZ = Path(config['working_dir']) / 'variant_filtering' / '04-dist' / 'vcf_gz-indexed.io'
    params:
        running_dir = Path(config['working_dir']) / 'variant_filtering' / '04-dist' / 'input'
    run:
        from camel.app.tools.bcftools.bcftoolsindex import BcftoolsIndex

        # Create working directory
        dir_working = Path(params.running_dir)
        dir_working.mkdir(parents=True, exist_ok=True)

        # Run tool
        bcftools_index = BcftoolsIndex(Camel.get_instance())
        step = Step(str(rule), bcftools_index, Camel.get_instance(), dir_working, config)
        SnakemakeUtils.add_pickle_inputs(bcftools_index, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_index, output)

rule variant_filtering_distance:
    """
    Filters the variants based on distance.
    """
    input:
        VCF_GZ = rules.variant_filtering_distance_index.output.VCF_GZ
    output:
        VCF_GZ = Path(config['working_dir']) / 'variant_filtering' / '04-dist' / 'vcf_gz.io',
        INFORMS = Path(config['working_dir']) / 'variant_filtering' / '04-dist' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'variant_filtering' / '04-dist',
        min_distance = variant_filtering.get_filtering_param(config, 'distance', 'min_distance'),
        keep_best = variant_filtering.get_filtering_param(config, 'distance', 'keep_best'),
        soft_filter = config['variant_filtering'].get('soft_filter', False)
    run:
        from camel.app.tools.variantfiltering.distancefilter import DistanceFilter
        distance_filter = DistanceFilter(Camel.get_instance())
        step = Step(str(rule), distance_filter, Camel.get_instance(), params.running_dir, config)
        if params.min_distance is not None:
            distance_filter.update_parameters(min_distance=params.min_distance)
        if params.keep_best is not None:
            distance_filter.update_parameters(keep_best=params.keep_best)
        distance_filter.update_parameters(soft_filter=params.soft_filter, seed=0)
        SnakemakeUtils.add_pickle_inputs(distance_filter, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(distance_filter, output)

rule variant_filtering_zscore_index:
    """
    Indexes the VCF file before applying the Z-score filter.
    """
    input:
        VCF_GZ = rules.variant_filtering_distance.output.VCF_GZ
    output:
        VCF_GZ = Path(config['working_dir']) / 'variant_filtering' / '05-zscore' / 'vcf_gz-indexed.io'
    params:
        running_dir = Path(config['working_dir']) / 'variant_filtering' / '05-zscore'
    run:
        from camel.app.tools.bcftools.bcftoolsindex import BcftoolsIndex
        bcftools_index  = BcftoolsIndex(Camel.get_instance())
        step = Step(str(rule), bcftools_index, Camel.get_instance(), params.running_dir, config)
        SnakemakeUtils.add_pickle_inputs(bcftools_index, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_index, output)

rule variant_filtering_zscore:
    """
    Filters variants based on the Z-score metric.
    """
    input:
        VCF_GZ = rules.variant_filtering_zscore_index.output.VCF_GZ,
        BAM = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_BAM if 'fasta_vcf' not in config['input_type'] else []
    output:
        VCF_GZ = Path(config['working_dir']) / 'variant_filtering' / '05-zscore' / 'vcf_gz.io',
        INFORMS = Path(config['working_dir']) / 'variant_filtering' / '05-zscore' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'variant_filtering' / '05-zscore',
        min_zscore = variant_filtering.get_filtering_param(config, 'zscore', 'min_zscore'),
        y_multiplier = variant_filtering.get_filtering_param(config, 'zscore', 'y_multiplier'),
        soft_filter = config['variant_filtering'].get('soft_filter', False)
    run:
        import logging
        import shutil
        from camel.app.tools.variantfiltering.zscorefilter import ZScoreFilter
        zscore_filter = ZScoreFilter(Camel.get_instance())
        step = Step(str(rule), zscore_filter, Camel.get_instance(), params.running_dir, config)
        if params.min_zscore is not None:
            zscore_filter.update_parameters(min_zscore=params.min_zscore)
        if params.y_multiplier is not None:
            zscore_filter.update_parameters(y_multiplier=params.y_multiplier)
        zscore_filter.update_parameters(soft_filter=params.soft_filter)
        if 'fasta_vcf' in config['input_type'] or len(SnakemakeUtils.load_object(Path(input.BAM))) == 0:
            logging.info("No BAM input found, skipping Z-score filter.")
            shutil.copyfile(input.VCF_GZ, output.VCF_GZ)
            SnakemakeUtils.dump_object({'variants_in': 'NA', 'variants_out': 'NA'}, Path(output.INFORMS))
        else:
            SnakemakeUtils.add_pickle_inputs(zscore_filter, input)
            step.run_step()
            SnakemakeUtils.dump_tool_outputs(zscore_filter, output)

rule variant_filtering_region:
    """
    Filters out regions that are problematic for variant calling based on an input BED file.
    """
    input:
        VCF_GZ = rules.variant_filtering_zscore.output.VCF_GZ
    output:
        VCF = Path(config['working_dir']) / 'variant_filtering' / '06-regions' / 'vcf.io',
        INFORMS = Path(config['working_dir']) / 'variant_filtering' / '06-regions' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'variant_filtering' / '06-regions',
        bed_file = variant_filtering.get_filtering_param(config, 'region', 'bed_file'),
        soft_filter = config['variant_filtering'].get('soft_filter', False)
    run:
        from camel.app.components.vcf.vcfutils import VCFUtils
        from camel.app.tools.bcftools.bcftoolsfilter import BcftoolsFilter
        from camel.app.io.tooliofile import ToolIOFile

        # Initialize tools
        bcftools_filter = BcftoolsFilter(Camel.get_instance())

        # Add input files
        vcf_file = SnakemakeUtils.load_object(Path(input.VCF_GZ))[0]
        bcftools_filter.add_input_files({'VCF_GZ': [vcf_file]})
        if params.bed_file is not None:
            bcftools_filter.add_input_files({'BED_exclude': [ToolIOFile(Path(params.bed_file))]})

        # Update parameters and run tool
        bcftools_filter.update_parameters(
            output_type='v', invert_targets=True, output_filename='variants_regions_removed.vcf',
            soft_filter=False if params.soft_filter is False else 'regions')
        step = Step(str(rule), bcftools_filter, Camel.get_instance(), params.running_dir, config)
        step.run_step()

        # Collect output
        SnakemakeUtils.dump_tool_output(bcftools_filter, 'VCF', Path(output.VCF))
        informs = {
            **bcftools_filter.informs,
            'variants_out': VCFUtils.count_variants(bcftools_filter.tool_outputs['VCF'][0].path),
            'variants_in': VCFUtils.count_variants(vcf_file.path),
        }
        SnakemakeUtils.dump_object(informs, Path(output.INFORMS))

rule variant_filtering_collect_stats:
    """
    Collects the stats for all of the filtering steps.
    """
    input:
        INFORMS_depth = rules.variant_filtering_depth.output.INFORMS,
        INFORMS_snp_qual = rules.variant_filtering_snp_quality.output.INFORMS,
        INFORMS_mapping_qual = rules.variant_filtering_mapping_quality.output.INFORMS,
        INFORMS_distance = rules.variant_filtering_distance.output.INFORMS,
        INFORMS_zscore = rules.variant_filtering_zscore.output.INFORMS,
        INFORMS_region = rules.variant_filtering_region.output.INFORMS
    output:
        JSON = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_STATS,
        INFORMS_ALL = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_INFORMS_ALL
    params:
        working_dir = Path(config['working_dir']) / 'variant_filtering' / 'stats'
    run:
        from camel.app.io.tooliofile import ToolIOFile
        import json

        filtering_data = {}
        all_informs = []
        for input_key in input.keys():
            filter_name = input_key.replace('INFORMS_', '')
            informs = SnakemakeUtils.load_object(Path(input[input_key]))
            filtering_data[filter_name] = informs
            if all([x in informs for x in ('_name', '_command')]):
                all_informs.append(informs)
        output_path = params.working_dir / 'stats.json'
        with output_path.open('w') as handle:
            json.dump(filtering_data, handle)
        SnakemakeUtils.dump_object([ToolIOFile(output_path)], Path(output.JSON))

        # Add tag to distinguish commands
        for inform in all_informs:
            inform['_tag'] = 'Variant filtering'
        SnakemakeUtils.dump_object(all_informs, Path(output.INFORMS_ALL))

rule variant_filtering_dump_summary_info:
    """
    Dumps the summary information for the variant filtering workflow.
    """
    input:
        JSON = rules.variant_filtering_collect_stats.output.JSON
    output:
        TSV = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_SUMMARY
    run:
        import json
        with open(SnakemakeUtils.load_object(Path(input.JSON))[0].path) as handle:
            filtering_info = json.load(handle)

        with open(output.TSV, 'w') as handle:
            for key, data in sorted(filtering_info.items()):
                handle.write('\t'.join([f'filt-{key}-in', str(data['variants_in'])]))
                handle.write('\n')
                handle.write('\t'.join([f'filt-{key}-out', str(data['variants_out'])]))
                handle.write('\n')
