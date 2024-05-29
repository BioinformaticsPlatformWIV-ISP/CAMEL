"""
This Snakefile performs the 51 SNP assay. It consists out of 3 assays:
- A positive control for M. tuberculosis complex members (SNP01)
- Genetic group based on katG codon 463 (SNP05), gyrA codon 95 (SNP06)
- Assign a best matching SNP cluster group (SCG) (SNP07-SNP51)
"""
from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import variant_calling, variant_filtering
from camel.scripts.mycobacteriumpipeline.snakefile import assay51snp


rule assay_51snp_init_db:
    """
    Converts the 51 SNP databases to CAMEL IO pickles. 
    """
    output:
        BED = Path(config['working_dir']) / '51snp' / 'bed.io',
        TSV = Path(config['working_dir']) / '51snp' / 'tsv.io'
    params:
        bed = config['51snp']['bed'],
        tsv = config['51snp']['profiles']
    run:
        from camel.app.io.tooliofile import ToolIOFile
        SnakemakeUtils.dump_object([ToolIOFile(Path(params.bed))], Path(output.BED))
        SnakemakeUtils.dump_object([ToolIOFile(Path(params.tsv))], Path(output.TSV))

rule assay_51snp_filter_vcf:
    """
    Extracts the SNPs at the 51 SNP positions.
    """
    input:
        VCF_GZ = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_VCF,
        BED = rules.assay_51snp_init_db.output.BED
    output:
        VCF = Path(config['working_dir']) / '51snp' / 'vcf.io'
    params:
        dir_ = Path(config['working_dir'], '51snp')
    run:
        from camel.app.tools.bcftools.bcftoolsfilter import BcftoolsFilter
        bcf_filter = BcftoolsFilter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(bcf_filter, input)
        bcf_filter.update_parameters(output_filename='51_snps.vcf')
        step = Step(str(rule), bcf_filter, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcf_filter, output)

rule assay_51snp_detect_info:
    """
    Detects the metadata for the 51SNP assay.
    """
    input:
        BED = rules.assay_51snp_init_db.output.BED,
        TSV = rules.assay_51snp_init_db.output.TSV,
        VCF = Path(config['working_dir']) / variant_calling.get_vcf(config),
        VCF_filt = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_VCF
    output:
        INFORMS = Path(config['working_dir']) / '51snp' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / '51snp'
    run:
        from camel.app.tools.pipelines.mycobacterium.assay51snpdetector import Assay51SnpDetector
        detector = Assay51SnpDetector(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(detector, input)
        step = Step(str(rule), detector, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(detector, output)

rule assay_51snp_report:
    """
    Creates the report for the 51 SNP assay.
    """
    input:
        INFORMS_detection = rules.assay_51snp_detect_info.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / assay51snp.OUTPUT_51SNP_REPORT
    params:
        dir_ = Path(config['working_dir']) / '51snp',
        sample_name = config['sample_name']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.mycobacterium.assay51snpreporter import Assay51SnpReporter
        spr = Assay51SnpReporter(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(spr, input)
        spr.add_input_files({'VAL_Sample': [ToolIOValue(params.sample_name)]})
        step = Step(str(rule), spr, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spr, output)

rule assay_51snp_report_empty:
    """
    Creates an empty report for the 51 SNP species identification.
    """
    output:
        VAL_HTML= Path(config['working_dir']) / assay51snp.OUTPUT_51SNP_REPORT_EMPTY
    params:
        dir_ = Path(config['working_dir']) / 'contamination_check' / 'report'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.pipelines.mycobacterium.assay51snpreporter import Assay51SnpReporter
        SnakePipelineUtils.create_empty_report_section(Assay51SnpReporter.TITLE, Path(output.VAL_HTML), 3)

rule assay_51snp_dump_summary_info:
    """
    Dumps the summary information from the 51SNP workflow.
    """
    input:
        INFORMS = rules.assay_51snp_detect_info.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / assay51snp.OUTPUT_51SNP_SUMMARY
    run:
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        summary_data = [
            ('51SNP-positive_control', informs['mtbc_pos_control']),
            ('51SNP-gyrB_group', informs['gyrB_group']),
            ('51SNP-genetic_group', informs['genetic_group']),
            ('51SNP-scg', informs['scg_profile'].scg),
            ('51SNP-st', informs['scg_profile'].st),
            ('51SNP-matching_snps', informs['scg_nb_snps_matched'])
        ]
        for i in range(1, 52):
            key = 'SNP{:02d}'.format(i)
            summary_data.append(('51SNP-{}'.format(key), informs[key]))

        with open(output.TSV, 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')
