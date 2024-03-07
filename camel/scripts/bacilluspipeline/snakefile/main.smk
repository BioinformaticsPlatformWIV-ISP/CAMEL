import shutil
from pathlib import Path

from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import core, assembly, downsampling, quast, confindr, trimming, trimming_illumina, \
    quality_checks, contamination_check_kraken, sequence_typing, amrfinder, trimming_ont, gene_detection, \
    mobsuite
from camel.scripts.bacilluspipeline.snakefile import btyper, ani

#######################
# Included Snakefiles #
#######################
include: core.SNAKEFILE_CORE
include: downsampling.SNAKEFILE_DOWNSAMPLING
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: trimming_ont.SNAKEFILE_TRIMMING_ONT
include: assembly.SNAKEFILE_ASSEMBLY
include: quast.SNAKEFILE_QUAST
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: confindr.SNAKEFILE_CONFINDR
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING
include: btyper.SNAKEFILE_BTYPER
include: amrfinder.SNAKEFILE_AMRFINDER
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: mobsuite.SNAKEFILE_MOB_SUITE
include: ani.SNAKEFILE_ANI

#########
# Rules #
#########
rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        HTML = config['output_report'],
        TSV = config['output_tabular']

#####################################
# Linking workflow inputs & outputs #
#####################################
rule link_fasta_to_tools_subtilis:
    """
    This rule links the output of the assembly workflow to the fastANI workflow if the species is B. subtilis.
    """
    input:
        FASTA = Path(config['working_dir'], assembly.OUTPUT_ASSEMBLY_FASTA)
    output:
        FASTA_ani = Path(config['working_dir']) / ani.INPUT_FASTA_ANI
    run:
        shutil.copyfile(Path(input.FASTA), Path(output.FASTA_ani))

rule link_fasta_to_tools_cereus:
    """
    This rule links the output of the assembly workflow to the BTyper workflow if the species is B. cereus.
    """
    input:
        FASTA = Path(config['working_dir'], assembly.OUTPUT_ASSEMBLY_FASTA)
    output:
        FASTA_btyper = Path(config['working_dir']) / btyper.INPUT_BTYPER_FASTA
    run:
        shutil.copyfile(Path(input.FASTA), Path(output.FASTA_btyper))

##########
# Report #
##########
rule report_init:
    """
    Creates the header section of the report.
    """
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-init.io'
    params:
        sample_name = config['sample_name'],
        config_input = config['input'],
        species = config['species'],
        input_type = config['input_type'],
        pipeline_info = config['pipeline'],
        detection_method = config['detection_method'],
        citation_keys = config['citations']
    run:
        import datetime
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.io.tooliovalue import ToolIOValue

        # Create the header section
        section = SnakePipelineUtils.create_input_section(
            sample_name=params.sample_name,
            date=datetime.datetime.now(),
            pipeline_version=params.pipeline_info['version'],
            input_files=ReportPipeline.format_input_string(params.config_input),
            input_type=params.input_type,
            extra_data=[('Detection method', params.detection_method),
             ('Species', f'<i>{params.species}</i>')],
            key_citation=params.citation_keys['main']
        )
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule report_create_commands_section:
    """
    Creates the section with the commands.
    """
    input:
        INFORMS_downsampling = downsampling.get_command_informs(config),
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_assembly = assembly.get_command_informs(config),
        INFORMS_quast = Path(config['working_dir']) / quast.OUTPUT_QUAST_INFORMS,
        INFORMS_busco = Path(config['working_dir']) / quast.OUTPUT_BUSCO_INFORMS,
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        # INFORMS_confindr = Path(config['working_dir']) / confindr.OUTPUT_CONFINDR_INFORMS if 'confindr' in config['analyses'] else [],
        INFORMS_assembly_map = assembly.get_qc_informs(config,config['input_type']),
        INFORMS_btyper = Path(config['working_dir']) / btyper.OUTPUT_INFORMS_BTYPER if 'btyper' in config['analyses'] else [],
        INFORMS_fastani = Path(config['working_dir']) / ani.OUTPUT_INFORMS_ANI if 'fastani' in config['analyses'] else [],
        INFORMS_amrfinder = Path(config['working_dir']) / str(amrfinder.OUTPUT_AMRFINDER_INFORMS) if 'amrfinder' in config['analyses'] else [],
        INFORMS_vfdb_core = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        INFORMS_plasmidfinder= Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else [],
        INFORMS_mob_suite = Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_INFORMS if 'mobsuite' in config['analyses'] else [],
        INFORMS_mlst_cereus = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='mlst_cereus') if 'mlst_cereus' in config['analyses'] else [],
        INFORMS_mlst_subtilis = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='mlst_subtilis') if 'mlst_subtilis' in config['analyses'] else []
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands-cereus.io'
    params:
        working_dir = config['working_dir']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        informs = []
        for content in [SnakemakeUtils.load_object(Path(io)) for io in input]:
            if type(content) is dict:
                informs.append(content)
            elif type(content) is list:
                informs.extend(content)
        section = SnakePipelineUtils.create_commands_section(informs, params.working_dir)
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule report_content_cereus:
    """
    Creates the main content of the report when the detected species is Bacillus cereus. 
    """
    input:
        report_init = rules.report_init.output.HTML,
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        report_quast = Path(config['working_dir']) / quast.OUTPUT_QUAST_REPORT,
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_confindr = confindr.get_report(config),
        report_adv_qc = Path(config['working_dir']) / str(quality_checks.OUTPUT_QUALITY_CHECKS_REPORT).format(
            input_type=config['input_type']),
        report_btyper = Path(config['working_dir']) / (btyper.OUTPUT_BTYPER_REPORT if 'btyper' in config['analyses'] else btyper.OUTPUT_BTYPER_REPORT_EMPTY),
        report_amrfinder = Path(config['working_dir']) / (amrfinder.OUTPUT_AMRFINDER_REPORT if 'amrfinder' in config['analyses'] else amrfinder.OUTPUT_AMRFINDER_REPORT_EMPTY),
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        report_plasmidfinder = gene_detection.get_gene_detection_report('plasmidfinder', config),
        report_mob_suite = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_REPORT if 'mobsuite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_REPORT_EMPTY),
        report_genomic_context = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_CONTEXT_REPORT if 'mobsuite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_CONTEXT_REPORT_EMPTY),
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        report_mlst = sequence_typing.get_sequence_typing_report('mlst_cereus', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst_cereus', config),
        report_commands = rules.report_create_commands_section.output.HTML,
        report_citations = Path(config['working_dir'],core.OUTPUT_HTML_CITATIONS)
    output:
        HTML = Path(config['working_dir']) / 'report' / 'report_cereus.html'
    params:
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        input_type = config['input_type']
    run:
        # Initialize the report
        Path(params.output_dir).mkdir(exist_ok=True, parents=True)
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakemakeUtils.load_object(Path(input.report_init))[0].value)
        report_structure = []

        # Core sections (shared)
        ReportPipeline.add_content_trim_basic_qc(
            report_structure, params.input_type, input.reports_downsampling, input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        ReportPipeline.add_content_contamination_check(
            report_structure, params.input_type, input.reports_contamination, input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))

        # Custom assays (B. cereus)
        report_structure.extend([
            ('BTyper3', 'btyper3', [Path(input.report_btyper)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (input.report_vfdb_core,)]),
            ('AMRFinder results', 'amrfinder', [Path(input.report_amrfinder)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_plasmidfinder, input.report_mob_suite, input.report_genomic_context)]),
            ('Sequence typing', 'st', [Path(x) for x in (input.report_rmlst, input.report_mlst, input.report_cgmlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ])
        SnakePipelineUtils.add_report_content(report, report_structure)
        report.save()

rule report_content_subtilis:
    """
    Creates the main content of the report when the detected species is Bacillus subtilis. 
    """
    input:
        report_init = rules.report_init.output.HTML,
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        report_quast = Path(config['working_dir']) / quast.OUTPUT_QUAST_REPORT,
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_confindr = confindr.get_report(config),
        report_adv_qc = Path(config['working_dir']) / str(quality_checks.OUTPUT_QUALITY_CHECKS_REPORT).format(
            input_type=config['input_type']),
        report_fastani = Path(config['working_dir']) / (ani.OUTPUT_ANI_REPORT if 'fastani' in config['analyses'] else ani.OUTPUT_ANI_REPORT_EMPTY),
        report_amrfinder = Path(config['working_dir']) / (amrfinder.OUTPUT_AMRFINDER_REPORT if 'amrfinder' in config['analyses'] else amrfinder.OUTPUT_AMRFINDER_REPORT_EMPTY),
        report_gmo = gene_detection.get_gene_detection_report('gmo', config),
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        report_plasmidfinder = gene_detection.get_gene_detection_report('plasmidfinder', config),
        report_mob_suite = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_REPORT if 'mobsuite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_REPORT_EMPTY),
        report_genomic_context = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_CONTEXT_REPORT if 'mobsuite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_CONTEXT_REPORT_EMPTY),
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        report_mlst = sequence_typing.get_sequence_typing_report('mlst_subtilis', config),
        report_citations = Path(config['working_dir'],core.OUTPUT_HTML_CITATIONS),
        report_commands = rules.report_create_commands_section.output.HTML
    output:
        HTML = Path(config['working_dir']) / 'report' / 'report_subtilis.html'
    params:
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        input_type = config['input_type']
    run:
        # Init report
        Path(params.output_dir).mkdir(exist_ok=True, parents=True)
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakemakeUtils.load_object(Path(input.report_init))[0].value)
        report_structure = []

        # Core sections
        ReportPipeline.add_content_trim_basic_qc(
            report_structure, params.input_type, input.reports_downsampling, input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        ReportPipeline.add_content_contamination_check(
            report_structure, params.input_type, input.reports_contamination, input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))

        # B. subtilis assays
        report_structure.extend([
            ('FastANI', 'fastani', [Path(input.report_fastani)]),
            ('GMO detection', 'gmo', [Path(input.report_gmo)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (input.report_vfdb_core,)]),
            ('AMRFinder results', 'amrfinder', [Path(input.report_amrfinder)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_plasmidfinder, input.report_mob_suite, input.report_genomic_context)]),
            ('Sequence typing', 'st', [Path(x) for x in (input.report_rmlst, input.report_mlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ])
        SnakePipelineUtils.add_report_content(report, report_structure)
        report.save()

rule report_select_by_species:
    """
    Selects the report content based on the detected species.
    """
    input:
        HTML = Path(config['working_dir']) / 'report' / f"report_{config['species']}.html"
    output:
        HTML = config['output_report']
    params:
        output_dir = config['output_dir']
    shell:
        """
        cp {input.HTML} {output.HTML};
        """

rule summary_combine_all:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        Path(config['working_dir'],core.OUTPUT_TSV_SUMMARY_INIT),
        downsampling.get_summaries(config),
        trimming.get_summaries(config),
        Path(config['working_dir']) / quast.OUTPUT_QUAST_SUMMARY,
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        contamination_check_kraken.get_summaries(config),
        confindr.get_summary(config),
        Path(config['working_dir']) / btyper.OUTPUT_BTYPER_SUMMARY if 'btyper' in config['analyses'] else [],
        Path(config['working_dir']) / amrfinder.OUTPUT_AMRFINDER_SUMMARY if 'amrfinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='gmo') if 'gmo' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else [],
        Path(config['working_dir']) / ani.OUTPUT_ANI_SUMMARY if 'fastani' in config['analyses'] else [],
        Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_SUMMARY if 'mobsuite' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='rmlst') if 'rmlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst_cereus') if 'mlst_cereus' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst_subtilis') if 'mlst_subtilis' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst_cereus') if 'cgmlst_cereus' in config['analyses'] else []
    output:
        TSV = config.get('output_tabular')
    run:
        with open(output.TSV, 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())

rule link_genomic_context:
    """
    Links the input databases to the genomic context assay.
    """
    input:
        # AMR
        TSV_amrfinder = Path(config['working_dir']) / 'amrfinder' / 'tsv.io' if 'amrfinder' in config['analyses'] else [],
        # Virulence
        TSV_gd_vfdb = Path(config['working_dir']) / 'gene_detection' / 'vfdb_core' / 'metadata' / 'tsv.io' if 'vfdb_core' in config['analyses'] else [],
        INFORMS_gd_vfdb = Path(config['working_dir']) / 'gene_detection' / 'vfdb_core' / 'db_manager' / 'informs.io' if 'vfdb_core' in config['analyses'] else []
    output:
        TSV = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'input' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'mob_suite' / 'genomic_context' / 'input' / 'informs.io'
    run:
        mobsuite.collect_genomic_context_input(input, Path(output.TSV), Path(output.INFORMS))
