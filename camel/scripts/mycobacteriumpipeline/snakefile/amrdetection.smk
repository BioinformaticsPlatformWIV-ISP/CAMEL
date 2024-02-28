import json
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import variant_calling, assembly, variant_filtering
from camel.scripts.mycobacteriumpipeline.snakefile import snplineage, amrdetection


rule amr_extract_variant_positions:
    """
    Extracts positions from the VCF file that are located in regions linked to AMR.
    """
    input:
        VCF_GZ = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ
    output:
        VCF = Path(config['working_dir']) / 'amr' / 'filtering' / 'vcf.io'
    params:
        dir_ = Path(config['working_dir']) / 'amr' / 'filtering',
        bed_regions = config['amr']['bed_regions']
    run:
        from camel.app.tools.bcftools.bcftoolsfilter import BcftoolsFilter
        bcf_filter = BcftoolsFilter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(bcf_filter, input)
        bcf_filter.add_input_files({'BED_include': [ToolIOFile(Path(params.bed_regions))]})
        step = Step(str(rule), bcf_filter, Camel.get_instance(), Path(params.dir_))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcf_filter, output)

rule amr_screen_mutations:
    """
    Screens the mutations detected in the AMR regions against the DB. 
    """
    input:
        VCF = rules.amr_extract_variant_positions.output.VCF,
        VCF_filt = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_VCF
    output:
        JSON = Path(config['working_dir']) / 'amr' / 'screen' / 'json.io',
        TSV = Path(config['working_dir']) / 'amr' / 'screen' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'amr' / 'screen' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'amr' / 'screen',
        db = config['amr']['mutation_db'],
        resistance_bed = config['amr']['bed_regions']
    run:
        from camel.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.pipelines.mycobacterium.amr.amrscreen import AMRScreen

        amr_screen = AMRScreen(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(amr_screen, input)
        amr_screen.add_input_files({
            'DB': [ToolIODirectory(Path(params.db))],
            'BED': [ToolIOFile(Path(params.resistance_bed))]
        })
        step = Step(str(rule), amr_screen, Camel.get_instance(), Path(params.dir_))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(amr_screen, output)

rule amr_export_positions:
    """
    Extracts the positions with potential AMR mutations and stores them in a text file (for determining base counts 
    with pileup).
    """
    input:
        VCF = rules.amr_extract_variant_positions.output.VCF
    output:
        TXT = Path(config['working_dir']) / 'amr' / 'filtering' / 'txt.io'
    params:
        dir_ = Path(config['working_dir']) / 'amr' / 'filtering'
    run:
        import vcf
        input_vcf = SnakemakeUtils.load_object(Path(input.VCF))[0].path
        output_path = Path(params.dir_, 'amr_positions.txt')
        with open(output_path, 'w') as handle_out, open(input_vcf) as handle_in:
            for variant in vcf.VCFReader(handle_in):
                handle_out.write('\t'.join([variant.CHROM, str(variant.POS)]))
                handle_out.write('\n')
        SnakemakeUtils.dump_object([ToolIOFile(output_path)], Path(output.TXT))

rule amr_pileup_variant_positions:
    """
    Creates a pileup for the variant positions. It is used to determine the ACTG counts.
    """
    input:
        FASTA = Path(config['working_dir']) / 'variant_calling' / 'reference' / 'fasta.io',
        BAM = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_BAM,
        TXT_POS = rules.amr_export_positions.output.TXT
    output:
        PILEUP = Path(config['working_dir']) / 'amr' / 'pileup' / 'pileup.io'
    params:
        dir_ = Path(config['working_dir']) / 'amr' / 'pileup'
    run:
        from camel.app.tools.samtools.samtoolsmpileup import SamtoolsMPileup
        samtools_mpileup = SamtoolsMPileup(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(samtools_mpileup, input)
        samtools_mpileup.update_parameters(count_orphans=True)
        step = Step(str(rule), samtools_mpileup, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_mpileup, output)

rule amr_parse_actg_counts:
    """
    Adds the ACTG counts to the annotated AMR mutations.
    """
    input:
        PILEUP = rules.amr_pileup_variant_positions.output.PILEUP
    output:
        JSON = Path(config['working_dir']) / 'amr' / 'pileup' / 'json_counts.io'
    run:
        from camel.app.components.mycobacterium import amrutils
        counts_by_pos = amrutils.parse_pileup(SnakemakeUtils.load_object(Path(input.PILEUP))[0].path)

        # Save output
        path_out = Path(output.JSON).parent / 'counts_by_pos.json'
        with path_out.open('w') as handle:
            json.dump(counts_by_pos, handle, indent=2)
        SnakemakeUtils.dump_object([ToolIOFile(path_out)], Path(output.JSON))

rule amr_predict_phenotype:
    """
    Predicts the phenotype for each antibiotic based on the detected mutations.
    """
    input:
        JSON = rules.amr_screen_mutations.output.JSON
    output:
        JSON = Path(config['working_dir']) / 'amr' / 'phenotype_prediction' / 'json_muts_by_ab.io'
    params:
        dir_ = Path(config['working_dir']) / 'amr' / 'phenotype_prediction',
        dir_amr_db = config['amr']['mutation_db']
    run:
        from camel.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.pipelines.mycobacterium.amr.amrphenotypepredictor import AMRPhenotypePredictor
        type_determination = AMRPhenotypePredictor(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(type_determination, input)
        type_determination.add_input_files({'DIR_DB': [ToolIODirectory(Path(params.dir_amr_db))]})
        step = Step(str(rule), type_determination, Camel.get_instance(), Path(params.dir_))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(type_determination, output)

rule amr_determine_resistance_type:
    """
    Determines the type of resistance (i.e. MDR, XDR) based on the predicted phenotypes.
    """
    input:
        JSON = rules.amr_predict_phenotype.output.JSON
    output:
        JSON = Path(config['working_dir']) / 'amr' / 'resistance_type' / 'json_amr_type.io'
    params:
        dir_ = Path(config['working_dir']) / 'amr' / 'resistance_type'
    run:
        from camel.app.tools.pipelines.mycobacterium.amr.amrtypedetermination import AMRTypeDetermination
        determination = AMRTypeDetermination(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(determination, input)
        step = Step(str(rule), determination, Camel.get_instance(), Path(params.dir_))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(determination, output)

rule amr_visualization_create_template:
    """
    Prepares the config file for the visualization of the mutations in AMR associated regions.
    """
    input:
        JSON = rules.amr_screen_mutations.output.JSON,
        TSV_depth = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_DEPTH_TSV
    output:
        TXT = Path(config['working_dir']) / 'amr' / 'visualization' / 'txt.io'
    params:
        dir_ = Path(config['working_dir']) / 'amr' / 'visualization',
        resistance_bed = config['amr']['bed_regions']
    run:
        from camel.app.tools.pipelines.mycobacterium.amr.amrcircostemplategeneration import AMRCircosTemplateGeneration
        templater = AMRCircosTemplateGeneration(Camel.get_instance())
        templater.update_parameters(spacing=500)
        templater.add_input_files({'BED': [ToolIOFile(Path(params.resistance_bed))]})
        SnakemakeUtils.add_pickle_inputs(templater, input)
        step = Step(str(rule), templater, Camel.get_instance(), Path(params.dir_))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(templater, output)

rule amr_visualization_circos:
    """
    Runs Circos to create the visualization.
    """
    input:
        TXT = rules.amr_visualization_create_template.output.TXT
    output:
        PNG = Path(config['working_dir']) / 'amr' / 'visualization' / 'png.io'
    params:
        dir_ = Path(config['working_dir']) / 'amr' / 'visualization'
    run:
        from camel.app.tools.circos.circos import Circos
        circos = Circos(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(circos, input)
        step = Step(str(rule), circos, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(circos, output)

rule amr_visualization_add_text:
    """
    Adds the text with the sample name, lineage and coverage information to the plot.
    """
    input:
        PNG = rules.amr_visualization_circos.output.PNG,
        INFORMS_coverage = Path(config['working_dir'] / assembly.get_depth_inform('fastq_pe')),
        INFORMS_lineage = Path(config['working_dir']) / snplineage.OUTPUT_SNP_LINEAGE_INFORMS
    output:
        PNG = Path(config['working_dir']) / 'amr' / 'visualization' / 'png-text.io'
    params:
        dir_ = Path(config['working_dir']) / 'amr' / 'visualization',
        sample_name = config['sample_name']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.mycobacterium.amr.amraddtext import AMRAddText
        text_adder = AMRAddText(Camel.get_instance())
        text_adder.add_input_files({'VAL_sample': [ToolIOValue(params.sample_name)]})
        SnakemakeUtils.add_pickle_inputs(text_adder, input)
        text_adder.run(Path(params.dir_))
        SnakemakeUtils.dump_tool_outputs(text_adder, output)

rule amr_create_report:
    """
    Creates a report for the variant based AMR detection.
    """
    input:
        JSON = rules.amr_screen_mutations.output.JSON,
        JSON_counts = rules.amr_parse_actg_counts.output.JSON,
        JSON_pheno = rules.amr_predict_phenotype.output.JSON,
        JSON_amr_type = rules.amr_determine_resistance_type.output.JSON,
        PNG = rules.amr_visualization_add_text.output.PNG,
        INFORMS_screen = rules.amr_screen_mutations.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / amrdetection.OUTPUT_AMR_REPORT
    params:
        dir_ = Path(config['working_dir']) / 'amr' / 'report',
        sample_name = config['sample_name'],
        bed_regions = config['amr']['bed_regions']
    run:
        from camel.app.tools.pipelines.mycobacterium.amr.amrreporter import AMRReporter
        reporter = AMRReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        reporter.add_input_files({'BED': [ToolIOFile(Path(params.bed_regions))]})
        step = Step(str(rule), reporter, Camel.get_instance(), Path(params.dir_))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule amr_create_report_empty:
    """
    Creates an empty report section for the AMR workflow when it is disabled.
    """
    output:
        HTML = Path(config['working_dir']) / amrdetection.OUTPUT_AMR_REPORT_EMPTY
    run:
        from camel.app.tools.pipelines.mycobacterium.amr.amrreporter import AMRReporter
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section(AMRReporter.TITLE, Path(output.HTML))

rule amr_dump_summary_info:
    """
    Dumps the summary information for the AMR assay.
    # TODO: Add DB version!
    """
    input:
        INFORMS_screening = rules.amr_screen_mutations.output.INFORMS,
        INFORMS_type = rules.amr_determine_resistance_type.output.JSON,
        INFORMS_pheno = rules.amr_predict_phenotype.output.JSON
    output:
        TSV = Path(config['working_dir']) / amrdetection.OUTPUT_AMR_SUMMARY
    run:
        import json
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        from camel.app.components.mycobacterium import amrutils

        output_data = []

        # DB version
        informs = SnakemakeUtils.load_object(Path(input.INFORMS_screening))
        output_data.append(['amr_db_version', informs['version']])

        # AMR type
        with open(SnakemakeUtils.load_object(Path(input.INFORMS_type))[0].path) as handle:
            data_pheno = json.load(handle)
        output_data.append(['amr_type', data_pheno['resistance_type']])
        output_data.append(['amr_first_line_resistant', data_pheno['first_line_resistant']])
        output_data.append(['amr_second_line_group_a_resistant', data_pheno['second_line_group_a_resistant']])
        output_data.append(['amr_second_line_group_b_resistant', data_pheno['second_line_group_b_resistant']])

        # Detected mutations
        with open(SnakemakeUtils.load_object(Path(input.INFORMS_pheno))[0].path) as handle:
            data_pheno = json.load(handle)
        for row in data_pheno:
            output_data.append([f"amr_pheno_{row['abbreviation']}", row['phenotype']])
            for level in amrutils.ConfidenceLevel:
                mutations = row['mutations'].get(level.value, [])
                output_data.append([
                    f"amr_mutations_{row['abbreviation']}_{level.value.replace(' ', '_')}",
                    ', '.join([f"{m['name']}" for m in mutations]) if len(mutations) > 0 else '-']
                )

        # Save output
        with open(output.TSV, 'w') as handle:
            for key, value in output_data:
                handle.write(f'{key}\t{str(value)}')
                handle.write('\n')
