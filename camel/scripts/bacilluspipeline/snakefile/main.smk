import shutil
from pathlib import Path

from camel.resources.snakefile import trimming, trimming_illumina, assembly_spades, \
    quality_checks, contamination_check_kraken, sequence_typing, downsampling, amrfinder, trimming_ont, gene_detection, \
    mobsuite, assembly_flye, medaka_polishing, short_read_polishing
from camel.scripts.bacilluspipeline.snakefile import btyper, ani

#######################
# Included Snakefiles #
#######################
# include: downsampling.SNAKEFILE_DOWNSAMPLING
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: trimming_ont.SNAKEFILE_TRIMMING_ONT
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: assembly_spades.SNAKEFILE_ASSEMBLY_SPADES
include: assembly_flye.SNAKEFILE_ASSEMBLY_FLYE
include: medaka_polishing.SNAKEFILE_MEDAKA_POLISHING
include: short_read_polishing.SNAKEFILE_POLISHING
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

####################################
# Linking workflow inputs & output #
####################################

rule prepare_fastq_inputs:
    """
    Prepares the fastq inputs
    """
    output:
        FASTQ_PE = Path(config['working_dir']) / 'illumina' / 'fastq.io' if config['read_type'] in ['illumina', 'hybrid'] else [],
        FASTQ_SE = Path(config['working_dir']) / 'nanopore' / 'fastq.io' if config['read_type'] in ['nanopore', 'hybrid'] else []
    params:
        config_input = config['input'],
        read_type = config['read_type']
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        if params.read_type == 'hybrid':
            SnakemakeUtils.dump_object([ToolIOFile(Path(x['path'])) for x in config['input']['fastq_pe']],
                Path(output.FASTQ_PE))
            SnakemakeUtils.dump_object([ToolIOFile(Path(x['path'])) for x in config['input']['fastq_se']],
                Path(output.FASTQ_SE))
        elif params.read_type == 'illumina':
            SnakemakeUtils.dump_object([ToolIOFile(Path(x['path'])) for x in config['input']['fastq_pe']],
                Path(output.FASTQ_PE))
        elif params.read_type == 'nanopore':
            SnakemakeUtils.dump_object([ToolIOFile(Path(x['path'])) for x in config['input']['fastq_se']],
                Path(output.FASTQ_SE))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

# rule link_downsampling_input_illumina:
#     """
#     Creates the FASTQ input for the downsampling step.
#     """
#     input:
#         Path(config['working_dir']) / 'input' / 'fastq_pe.io' if config['read_type'] in ['hybrid', 'illumina'] else []
#     output:
#         FASTQ = Path(config['working_dir']) / downsampling.INPUT_DOWNSAMPLING_FASTQ
#     params:
#         config_input = config['input'],
#         read_type = config['read_type']
#     run:
#
#         from camel.app.io.tooliofile import ToolIOFile
#         from camel.app.snakemake.snakemakeutils import SnakemakeUtils
#         if params.read_type in ['illumina', 'hybrid']:
#             SnakemakeUtils.dump_object([ToolIOFile(Path(x['path'])) for x in config['input']['fastq_pe']],
#                 Path(output.FASTQ))
#         elif params.read_type == 'nanopore':
#             SnakemakeUtils.dump_object([ToolIOFile(Path(x['path'])) for x in config['input']['fastq_se']],
#                 Path(output.FASTQ))
#         else:
#             raise ValueError(f'Unsupported read type: {params.read_type}')
#
# rule link_downsampling_input_nanopore:
#     """
#     Creates the FASTQ input for the downsampling step.
#     """
#     output:
#         FASTQ = Path(config['working_dir']) / downsampling.INPUT_DOWNSAMPLING_FASTQ
#     params:
#         config_input = config['input'],
#         read_type = config['read_type']
#     run:
#         from camel.app.io.tooliofile import ToolIOFile
#         from camel.app.snakemake.snakemakeutils import SnakemakeUtils
#         if params.read_type in ['illumina', 'hybrid']:
#             SnakemakeUtils.dump_object([ToolIOFile(Path(x['path'])) for x in config['input']['fastq_pe']],
#                 Path(output.FASTQ))
#         elif params.read_type in ['nanopore', 'hybrid']:
#             SnakemakeUtils.dump_object([ToolIOFile(Path(x['path'])) for x in config['input']['fastq_se']],
#                 Path(output.FASTQ))
#         else:
#             raise ValueError(f'Unsupported read type: {params.read_type}')
#
rule link_downsampling_to_trimming_workflows:
    """
    Links the downsampling output to the input of the trimming workflows.  
    """
    input:
        FASTQ_PE = rules.prepare_fastq_inputs.output.FASTQ_PE if config['read_type'] in ['illumina', 'hybrid'] else [],
        FASTQ_SE = rules.prepare_fastq_inputs.output.FASTQ_SE if config['read_type'] in ['nanopore', 'hybrid'] else []
        # FASTQ = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_FASTQ
    output:
        FASTQ_ilmn = Path(config['working_dir']) / trimming_illumina.INPUT_TRIMMOMATIC_FASTQ if config['read_type'] in ['illumina', 'hybrid'] else [],
        FASTQ_ont = Path(config['working_dir']) / trimming_ont.INPUT_ONT_FASTQ if config['read_type'] in ['nanopore', 'hybrid'] else []
    params:
        read_type = config['read_type']
    run:
        if params.read_type in ['nanopore', 'hybrid']:
            shutil.copyfile(Path(input.FASTQ_SE),Path(output.FASTQ_ont))
        if params.read_type in ['illumina', 'hybrid']:
            shutil.copyfile(Path(input.FASTQ_PE),Path(output.FASTQ_ilmn))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

rule link_fastq_to_fq_dict:
    """
    Creates an IO object with the trimmed FASTQ files.
    Other workflows such as Kraken or de novo assembly rely on this dictionary to get input files (PE or SE).
    """
    input:
        FASTQ_ilmn = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_DICT if config['read_type'] in ['illumina', 'hybrid'] else [],
        FASTQ_ont = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_DICT if config['read_type'] in ['nanopore', 'hybrid'] else []
    output:
        IO_FASTQ = Path(config['working_dir']) / 'fq_dict.io'
    params:
        read_type = config['read_type']
    run:
        if params.read_type == 'illumina':
            shutil.copyfile(Path(input.FASTQ_ilmn), Path(output.IO_FASTQ))
        elif params.read_type == 'nanopore':
            shutil.copyfile(Path(input.FASTQ_ont), Path(output.IO_FASTQ))
        elif params.read_type == 'hybrid':
            nanopore_reads = SnakemakeUtils.load_object(Path(input.FASTQ_ont))
            illumina_reads = SnakemakeUtils.load_object(Path(input.FASTQ_ilmn))
            nanopore_reads.update(illumina_reads)
            SnakemakeUtils.dump_object(nanopore_reads, Path(output.IO_FASTQ))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

rule copy_flye_to_medaka:
    """
    This rule copies necessary files to the short read polishing snakemake.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_flye.OUTPUT_ASSEMBLY_FASTA
    output:
        FASTA = Path(config['working_dir']) / str(medaka_polishing.INPUT_ASSEMBLY_FASTA).format(assembly_type=config.get('assembly_type', 'long_read_assembly'))
    run:
        shutil.copyfile(input.FASTA, output.FASTA)

rule copy_fasta_to_polishing:
    """
    This rule copies necessary files to the short read polishing snakemake.
    """
    input:
        FASTA = Path(config['working_dir']) / str(medaka_polishing.OUTPUT_ASSEMBLY_FASTA).format(assembly_type=config.get('assembly_type', 'long_read_assembly'))
    output:
        FASTA = Path(config['working_dir']) / str(short_read_polishing.INPUT_ASSEMBLY_FASTA).format(assembly_type=config.get('assembly_type', 'long_read_assembly'))
    run:
        shutil.copyfile(input.FASTA,output.FASTA)

rule select_fasta_file:
    """
    This rule selects the fasta file to send to other workflows.
    """
    input:
        FASTA_spades = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA if config['read_type'] == 'illumina' else [],
        FASTA_medaka = Path(config['working_dir']) / str(medaka_polishing.OUTPUT_ASSEMBLY_FASTA).format(assembly_type=config.get('assembly_type', 'long_read_assembly')) if config['read_type'] == 'nanopore' else [],
        FASTA_polished = Path(config['working_dir']) / str(short_read_polishing.OUTPUT_POLISHING_FASTA).format(assembly_type=config.get('assembly_type', 'long_read_assembly')) if config['read_type'] == 'hybrid' else []
    output:
        FASTA = Path(config['working_dir']) / 'fasta.io'
    params:
        read_type=config['read_type']
    run:
        if params.read_type == 'illumina':
            shutil.copyfile(Path(input.FASTA_spades),Path(output.FASTA))
        elif params.read_type == 'nanopore':
            shutil.copyfile(Path(input.FASTA_medaka),Path(output.FASTA))
        elif params.read_type == 'hybrid':
            shutil.copyfile(Path(input.FASTA_polished), output.FASTA)
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

rule link_fasta_to_gene_detection:
    """
    This rule links the output of the assembly workflows to the gene detection workflow.
    """
    input:
        FASTA = rules.select_fasta_file.output.FASTA
    output:
        FASTA_genedetection = Path(config['working_dir']) / gene_detection.INPUT_GENE_DETECTION_FASTA
    run:
        shutil.copyfile(Path(input.FASTA), Path(output.FASTA_genedetection))

rule link_fasta_to_typing:
    """
    This rule links the output of the assembly workflows to the sequence typing workflow.
    """
    input:
        FASTA = rules.select_fasta_file.output.FASTA
    output:
        FASTA_typing = Path(config['working_dir']) / sequence_typing.INPUT_FASTA
    run:
        shutil.copyfile(Path(input.FASTA), output.FASTA_typing)

rule link_fasta_to_tools_all:
    """
    This rule links the output of the assembly workflow to the amrfinder and mobsuite workflows.
    """
    input:
        FASTA = rules.select_fasta_file.output.FASTA
    output:
        FASTA_amrfinder = Path(config['working_dir']) / amrfinder.INPUT_AMRFINDER_FASTA,
        FASTA_mobsuite = Path(config['working_dir']) / mobsuite.INPUT_MOBSUITE_FASTA
    run:
        shutil.copyfile(Path(input.FASTA),Path(output.FASTA_amrfinder))
        shutil.copyfile(Path(input.FASTA),Path(output.FASTA_mobsuite))

rule link_fasta_to_tools_subtilis:
    """
    This rule links the output of the assembly workflow to the fastANI workflow if the species is B. subtilis.
    """
    input:
        FASTA = rules.select_fasta_file.output.FASTA
    output:
        FASTA_ani = Path(config['working_dir']) / ani.INPUT_FASTA_ANI
    run:
        shutil.copyfile(Path(input.FASTA),Path(output.FASTA_ani))

rule link_fasta_to_tools_cereus:
    """
    This rule links the output of the assembly workflow to the BTyper workflow if the species is B. cereus.
    """
    input:
        FASTA = rules.select_fasta_file.output.FASTA
    output:
        FASTA_btyper = Path(config['working_dir']) / btyper.INPUT_BTYPER_FASTA
    run:
        shutil.copyfile(Path(input.FASTA),Path(output.FASTA_btyper))

#############
# Read type #
#############

rule select_assembly_output:
    """
    Selects the assembly output based on the read type.
    """
    input:
        HTML = [Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_REPORT if config['read_type'] == 'illumina' else [],
                Path(config['working_dir']) / assembly_flye.OUTPUT_ASSEMBLY_REPORT if config['read_type'] == 'nanopore' else [],
                Path(config['working_dir']) / str(short_read_polishing.OUTPUT_ASSEMBLY_REPORT).format(assembly_type=config.get('assembly_type', 'long_read_assembly')) if config['read_type'] == 'hybrid' else []],
        TSV = [Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_SUMMARY if config['read_type'] == 'illumina' else [],
               Path(config['working_dir']) / assembly_flye.OUTPUT_ASSEMBLY_SUMMARY if config['read_type'] == 'nanopore' else [],
               Path(config['working_dir']) / str(short_read_polishing.OUTPUT_ASSEMBLY_SUMMARY).format(assembly_type=config.get('assembly_type', 'long_read_assembly')) if config['read_type'] == 'hybrid' else []],
        INFORMS = [Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_INFORMS if config['read_type'] == 'illumina' else [],
                   Path(config['working_dir']) / assembly_flye.OUTPUT_ASSEMBLY_INFORMS if config['read_type'] == 'nanopore' else [],
                   Path(config['working_dir']) / str(short_read_polishing.OUTPUT_ASSEMBLY_INFORMS).format(assembly_type=config.get('assembly_type', 'long_read_assembly')) if config['read_type'] == 'hybrid' else []],
        INFORMS_filt = [Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FILTERING_INFORMS if config['read_type'] == 'illumina' else [],
                        Path(config['working_dir']) / assembly_flye.OUTPUT_ASSEMBLY_FILTERING_INFORMS if config['read_type'] == 'nanopore' else [],
                        Path(config['working_dir']) / str(short_read_polishing.OUTPUT_ASSEMBLY_FILTERING_INFORMS).format(assembly_type=config.get('assembly_type', 'long_read_assembly')) if config['read_type'] == 'hybrid' else []]
    output:
        HTML = Path(config['working_dir']) / 'read_type' / 'assembly' / 'html.io',
        TSV = Path(config['working_dir']) / 'read_type' / 'assembly' / 'summary.tsv',
        INFORMS = Path(config['working_dir']) / 'read_type' / 'assembly' / 'informs.io',
        INFORMS_filt = Path(config['working_dir']) / 'read_type' / 'assembly' / 'informs-filt.io'
    run:
        shutil.copyfile([f for f in input.HTML if f][0], Path(output.HTML))
        shutil.copyfile([f for f in input.TSV if f][0], Path(output.TSV))
        shutil.copyfile([f for f in input.INFORMS if f][0], Path(output.INFORMS))
        shutil.copyfile([f for f in input.INFORMS_filt if f][0], Path(output.INFORMS_filt))

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
        read_type = config['read_type'],
        pipeline_info = config['pipeline'],
        detection_method = config['detection_method'],
        citation_keys = config['citations']
    run:
        import datetime
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.io.tooliovalue import ToolIOValue

        # Create header section
        section = SnakePipelineUtils.create_input_section(
            params.sample_name,
            datetime.datetime.now(),
            params.pipeline_info['version'], ', '.join(
                input_file['name'] for _, input_files in params.config_input.items() for input_file in input_files),
            [('Detection method', params.detection_method),
             ('Read type', str(params.read_type)),
             ('Species', f'<i>{params.species}</i>')],
            params.citation_keys['main']
        )
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule report_pickle_citations:
    """
    This rule creates a pickle with a report section containing the citations.
    """
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-citations.io'
    params:
        citation_keys = config['citations']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        section = SnakePipelineUtils.create_citations_section(
            params.citation_keys['other'], params.citation_keys['main'])
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule report_create_commands_section:
    """
    Creates the section with the commands.
    """
    input:
        # INFORMS_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_INFORMS,
        INFORMS_trimming_illumina = trimming.get_trimming_command_informs(config, 'illumina') if config['read_type'] in ['illumina', 'hybrid'] else [],
        INFORMS_trimming_nanopore = trimming.get_trimming_command_informs(config, 'nanopore') if config['read_type'] in ['nanopore', 'hybrid'] else [],
        INFORMS_assembly = rules.select_assembly_output.output.INFORMS,
        INFORMS_assembly_filt = rules.select_assembly_output.output.INFORMS_filt,
        INFORMS_kraken_illumina = contamination_check_kraken.get_contamination_check_kraken_informs(config, 'illumina') if config['read_type'] in ['illumina', 'hybrid'] else [],
        INFORMS_kraken_nanopore = contamination_check_kraken.get_contamination_check_kraken_informs(config, 'nanopore') if config['read_type'] in ['nanopore', 'hybrid'] else [],
        INFORMS_mapping_illumina = quality_checks.get_mapping_rate_informs(config, 'illumina') if config['read_type'] in ['illumina', 'hybrid'] else [],
        INFORMS_mapping_nanopore = quality_checks.get_mapping_rate_informs(config, 'nanopore') if config['read_type'] in ['nanopore', 'hybrid'] else [],
        INFORMS_depth_illumina = quality_checks.get_depth_informs(config, 'illumina') if config['read_type'] in ['illumina', 'hybrid'] else [],
        INFORMS_depth_nanopore = quality_checks.get_depth_informs(config, 'nanopore') if config['read_type'] in ['nanopore', 'hybrid'] else [],
        INFORMS_btyper = Path(config['working_dir']) / btyper.OUTPUT_INFORMS_BTYPER if 'btyper' in config['analyses'] else [],
        INFORMS_fastani = Path(config['working_dir']) / ani.OUTPUT_INFORMS_ANI if 'fastani' in config['analyses'] else [],
        INFORMS_amrfinder = Path(config['working_dir']) / str(amrfinder.OUTPUT_AMRFINDER_INFORMS) if 'amrfinder' in config['analyses'] else [],
        INFORMS_vfdb_core = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        INFORMS_plasmidfinder= Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else [],
        INFORMS_mob_suite = Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_INFORMS if 'mob_suite' in config['analyses'] else [],
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
        # report_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_REPORT,
        report_trimming_illumina = trimming.get_trimming_report(config, 'illumina') if config['read_type'] in ['illumina', 'hybrid'] else [],
        report_trimming_nanopore = trimming.get_trimming_report(config,'nanopore') if config['read_type'] in ['nanopore', 'hybrid'] else [],
        report_assembly = rules.select_assembly_output.output.HTML,
        report_kraken_illumina = contamination_check_kraken.get_contamination_check_report(config, 'illumina') if config['read_type'] in ['illumina', 'hybrid'] else [],
        report_kraken_nanopore = contamination_check_kraken.get_contamination_check_report(config, 'nanopore') if config['read_type'] in ['nanopore', 'hybrid'] else [],
        report_adv_qc_illumina = quality_checks.get_qc_report(config, 'illumina') if config['read_type'] in ['illumina','hybrid'] else [],
        report_adv_qc_nanopore = quality_checks.get_qc_report(config, 'nanopore') if config['read_type'] in ['nanopore','hybrid'] else [],
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
        report_citations = rules.report_pickle_citations.output.HTML
    output:
        HTML = Path(config['working_dir']) / 'report' / 'report_cereus.html'
    params:
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline']
    run:
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        Path(params.output_dir).mkdir(exist_ok=True)
        report.add_html_object(SnakemakeUtils.load_object(Path(input.report_init))[0].value)
        report_structure = [
            ('Read trimming and basic QC', 'trim', [Path(x) for x in (input.report_trimming_illumina, input.report_trimming_nanopore) if x != []]),
            ('Assembly', 'assembly', [Path(input.report_assembly)]),
            ('Advanced QC', 'adv_qc', [Path(x) for x in (input.report_adv_qc_illumina, input.report_adv_qc_nanopore, input.report_kraken_illumina, input.report_kraken_nanopore) if x != []]),
            ('BTyper3', 'btyper3', [Path(input.report_btyper)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (input.report_vfdb_core,)]),
            ('AMRFinder results', 'amrfinder', [Path(input.report_amrfinder)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_plasmidfinder, input.report_mob_suite, input.report_genomic_context)]),
            ('Sequence typing', 'st', [Path(x) for x in (input.report_rmlst, input.report_mlst, input.report_cgmlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ]
        SnakePipelineUtils.add_report_content(report, report_structure)
        report.save()

rule report_content_subtilis:
    """
    Creates the main content of the report when the detected species is Bacillus subtilis. 
    """
    input:
        report_init = rules.report_init.output.HTML,
        # report_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_REPORT,
        report_trimming_illumina = trimming.get_trimming_report(config,'illumina') if config['read_type'] in['illumina','hybrid'] else [],
        report_trimming_nanopore = trimming.get_trimming_report(config,'nanopore') if config['read_type'] in['nanopore','hybrid'] else [],
        report_assembly = rules.select_assembly_output.output.HTML,
        report_kraken_illumina = contamination_check_kraken.get_contamination_check_report(config, 'illumina') if (config['read_type'] in ['illumina', 'hybrid'] and 'kraken' in config['analyses']) else [],
        report_kraken_nanopore = contamination_check_kraken.get_contamination_check_report(config, 'nanopore') if (config['read_type'] in ['nanopore', 'hybrid'] and 'kraken' in config['analyses']) else [],
        report_adv_qc_illumina = quality_checks.get_qc_report(config, 'illumina') if config['read_type'] in ['illumina', 'hybrid'] else [],
        report_adv_qc_nanopore = quality_checks.get_qc_report(config, 'nanopore') if config['read_type'] in ['nanopore', 'hybrid'] else [],
        report_fastani = Path(config['working_dir']) / (ani.OUTPUT_ANI_REPORT if 'fastani' in config['analyses'] else ani.OUTPUT_ANI_REPORT_EMPTY),
        report_amrfinder = Path(config['working_dir']) / (amrfinder.OUTPUT_AMRFINDER_REPORT if 'amrfinder' in config['analyses'] else amrfinder.OUTPUT_AMRFINDER_REPORT_EMPTY),
        report_gmo = gene_detection.get_gene_detection_report('gmo', config),
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        report_plasmidfinder = gene_detection.get_gene_detection_report('plasmidfinder', config),
        report_mob_suite = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_REPORT if 'mobsuite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_REPORT_EMPTY),
        report_genomic_context = Path(config['working_dir']) / (mobsuite.OUTPUT_MOB_SUITE_CONTEXT_REPORT if 'mobsuite' in config['analyses'] else mobsuite.OUTPUT_MOB_SUITE_CONTEXT_REPORT_EMPTY),
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        report_mlst = sequence_typing.get_sequence_typing_report('mlst_subtilis', config),
        report_citations = rules.report_pickle_citations.output.HTML,
        report_commands = rules.report_create_commands_section.output.HTML
    output:
        HTML = Path(config['working_dir']) / 'report' / 'report_subtilis.html'
    params:
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline']
    run:
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakemakeUtils.load_object(Path(input.report_init))[0].value)
        report_structure = [
            # ('Read trimming and basic QC', 'trim', [Path(input.report_downsampling), Path(input.report_trimming)]),
            ('Read trimming and basic QC', 'trim', [Path(x) for x in (input.report_trimming_illumina, input.report_trimming_nanopore) if x != []]),
            ('Assembly', 'assembly', [Path(input.report_assembly)]),
            ('Advanced QC', 'adv_qc', [Path(x) for x in (input.report_adv_qc_illumina, input.report_adv_qc_nanopore, input.report_kraken_illumina, input.report_kraken_nanopore) if x != []]),
            ('FastANI', 'fastani', [Path(input.report_fastani)]),
            ('GMO detection', 'gmo', [Path(input.report_gmo)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (input.report_vfdb_core,)]),
            ('AMRFinder results', 'amrfinder', [Path(input.report_amrfinder)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_plasmidfinder, input.report_mob_suite, input.report_genomic_context)]),
            ('Sequence typing', 'st', [Path(x) for x in (input.report_rmlst, input.report_mlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ]
        SnakePipelineUtils.add_report_content(report,report_structure)
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
        cp {input.HTML} {output.HTML}
        """

rule summary_init:
    """
    Initializes the summary output file.
    """
    output:
        TSV = Path(config['working_dir']) / 'summary' / 'summary-init.tsv'
    run:
        import datetime
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        analysis_date = datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT)
        input_filenames = ', '.join(
            input_file['name'] for _, input_files in config['input'].items() for input_file in input_files)
        with open(output.TSV, 'w') as handle:
            for kv_pair in [
                ('pipeline_name', config['pipeline']['name']),
                ('pipeline_version', config['pipeline']['version']),
                ('sample', config['sample_name']),
                ('input_files', input_filenames),
                ('analysis_date', analysis_date),
                ('detection_method', config['detection_method'])]:
                handle.write('\t'.join(kv_pair))
                handle.write('\n')

rule summary_combine_all:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        rules.summary_init.output.TSV,
        # Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_SUMMARY,
        trimming.get_trimming_summary(config, 'illumina') if config['read_type'] in ['illumina', 'hybrid'] else [],
        trimming.get_trimming_summary(config, 'nanopore') if config['read_type'] in ['nanopore', 'hybrid'] else [],
        Path(config['working_dir']) / rules.select_assembly_output.output.TSV,
        Path(config['working_dir']) / medaka_polishing.OUTPUT_ASSEMBLY_SUMMARY if config['read_type'] == 'nanopore' else [],
        Path(config['working_dir']) / str(quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY).format(read_type='illumina') if config['read_type'] in ['illumina', 'hybrid'] else [],
        Path(config['working_dir']) / str(quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY).format(read_type='nanopore') if config['read_type'] in ['nanopore', 'hybrid'] else [],
        Path(config['working_dir']) / str(contamination_check_kraken.OUTPUT_CONTAMINATION_SUMMARY).format(read_type='illumina') if (config['read_type'] in ['illumina', 'hybrid'] and 'kraken' in config['analyses']) else [],
        Path(config['working_dir']) / str(contamination_check_kraken.OUTPUT_CONTAMINATION_SUMMARY).format(read_type='nanopore') if (config['read_type'] in ['nanopore', 'hybrid'] and 'kraken' in config['analyses']) else [],
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
        config.get('output_tabular')
    run:
        with open(output[0], 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
