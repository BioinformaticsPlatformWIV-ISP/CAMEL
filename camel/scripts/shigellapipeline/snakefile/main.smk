from pathlib import Path

from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import trimming, trimming_illumina, assembly_spades, \
    quality_checks, contamination_check_kraken, gene_detection, pointfinder, variant_calling, variant_filtering, \
    sequence_typing, downsampling
from camel.scripts.shigellapipeline.snakefile import subspecies_identification, flexneritype

#######################
# Included Snakefiles #
#######################
include: downsampling.SNAKEFILE_DOWNSAMPLING
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: assembly_spades.SNAKEFILE_ASSEMBLY_SPADES
include: variant_calling.SNAKEFILE_VARIANT_CALLING
include: variant_filtering.SNAKEFILE_VARIANT_FILTERING
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: pointfinder.SNAKEFILE_POINTFINDER
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING
include: subspecies_identification.SNAKEFILE_SUBSPECIES
include: flexneritype.SNAKEFILE_FLEXNERITYPE

#########
# Rules #
#########
rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        config['output_report'],
        config['output_tabular']

rule link_downsampling_input:
    """
    Creates the FASTQ input for the downsampling step. 
    """
    input:
        FASTQ_PE = [entry['path'] for entry in config['input']['fastq_pe']]
    output:
        FASTQ = Path(config['working_dir']) / downsampling.INPUT_DOWNSAMPLING_FASTQ
    run:
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        from camel.app.io.tooliofile import ToolIOFile
        SnakemakeUtils.dump_object([ToolIOFile(Path(x)) for x in input.FASTQ_PE], Path(output.FASTQ))

rule link_trimmomatic_input:
    """
    Links the downsampling output to the input of the trimmomatic workflow.  
    """
    input:
        FASTQ = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_FASTQ
    output:
        FASTQ = Path(config['working_dir']) / trimming_illumina.INPUT_TRIMMOMATIC_FASTQ
    shell:
        """
        cp {input.FASTQ} {output.FASTQ};
        """

rule init_summary:
    """
    Initializes the summary output file.
    """
    output:
        TSV = Path(config['working_dir']) / 'summary' / 'summary-init.tsv'
    run:
        import datetime
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

rule select_fastq:
    """
    This rule creates an IO object with the trimmed FASTQ files.
    Other workflows such as Kraken or Assembly rely on this dictionary to get input files (PE or SE).
    """
    input:
        FASTQ_PE = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_DICT
    output:
        IO_FASTQ = Path(config['working_dir']) / 'fq_dict.io'
    shell:
        "cp {input.FASTQ_PE} {output.IO_FASTQ};"

rule select_fasta:
    """
    This rules links the output of the assembly workflow to the other workflows. 
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA
    output:
        FASTA = Path(config['working_dir']) / gene_detection.INPUT_GENE_DETECTION_FASTA
    shell:
        "cp {input.FASTA} {output.FASTA};"

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

rule report_command_section:
    input:
        INFORMS_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_INFORMS,
        INFORMS_trimming = trimming.get_trimming_command_informs(config),
        INFORMS_assembly = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_INFORMS,
        INFORMS_assembly_filt=Path(config['working_dir']) / 'assembly_spades' / 'filtering' / 'informs.io',
        INFORMS_kraken = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS if 'kraken' in config['analyses'] else [],
        INFORMS_mapping = quality_checks.get_mapping_rate_informs(config),
        INFORMS_depth = quality_checks.get_depth_informs(config),
        INFORMS_variant_calling_all = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_INFORMS_ALL,
        INFORMS_variant_filtering_all = Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_INFORMS_ALL,
        INFORMS_resfinder = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='resfinder') if 'resfinder' in config['analyses'] else [],
        INFORMS_argannot = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='argannot') if 'argannot' in config['analyses'] else [],
        INFORMS_card = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='card') if 'card' in config['analyses'] else [],
        INFORMS_ncbi_amr = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
        INFORMS_pointfinder = Path(config['working_dir']) / pointfinder.OUTPUT_POINTFINDER_INFORMS if 'pointfinder' in config['analyses'] else [],
        INFORMS_virulence = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        INFORMS_virulence_shiga = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='virulencefinder_shiga') if 'virulencefinder' in config['analyses'] else [],
        INFORMS_plasmidfinder = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else []
    output:
        VAL_HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        working_dir=config['working_dir']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        informs = []
        for content in [SnakemakeUtils.load_object(Path(io)) for io in input]:
            if type(content) is dict:
                informs.append(content)
            elif type(content) is list:
                informs.extend(content)
        section = SnakePipelineUtils.create_commands_section(informs, Path(params.working_dir))
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.VAL_HTML))

rule combine_reports:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        report_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_REPORT,
        report_trimming = trimming.get_trimming_report(config),
        report_assembly = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_REPORT,
        report_kraken = Path(config['working_dir']) / (contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT if 'kraken' in config['analyses'] else contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY),
        report_adv_qc = Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT,
        report_variant = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_REPORT,
        report_pointfinder = Path(config['working_dir']) / (pointfinder.OUTPUT_POINTFINDER_REPORT if 'pointfinder' in config['analyses'] else pointfinder.OUTPUT_POINTFINDER_REPORT_EMPTY),
        # Gene detection
        report_resfinder = gene_detection.get_gene_detection_report('resfinder', config),
        report_argannot = gene_detection.get_gene_detection_report('argannot', config),
        report_card = gene_detection.get_gene_detection_report('card', config),
        report_ncbi_amr = gene_detection.get_gene_detection_report('ncbi_amr', config),
        report_virulence = gene_detection.get_gene_detection_report('virulencefinder', config),
        report_virulence_shiga = gene_detection.get_gene_detection_report('virulencefinder_shiga', config, 'virulencefinder'),
        report_plasmidfinder = gene_detection.get_gene_detection_report('plasmidfinder', config),
        # Shigella typing
        report_species = Path(config['working_dir']) / (subspecies_identification.OUTPUT_SPECIES_REPORT if 'identification' in config['analyses'] else subspecies_identification.OUTPUT_SPECIES_REPORT_EMPTY),
        report_subspecies = Path(config['working_dir']) /  (subspecies_identification.OUTPUT_SUBSPECIES_REPORT if 'identification' in config['analyses'] else subspecies_identification.OUTPUT_SUBSPECIES_REPORT_EMPTY),
        report_flexneri = Path(config['working_dir']) / (flexneritype.OUTPUT_FLEXNERI_REPORT if 'identification' in config['analyses'] else flexneritype.OUTPUT_FLEXNERI_REPORT_EMPTY),
        report_subspecies_db = Path(config['working_dir']) / 'subspecies_identification' / 'report' / 'html-db.io',
        # Sequence typing
        report_mlst_warwick = sequence_typing.get_sequence_typing_report('mlst_warwick', config),
        report_mlst_pasteur = sequence_typing.get_sequence_typing_report('mlst_pasteur', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        report_citations = rules.report_pickle_citations.output.HTML,
        report_commands = rules.report_command_section.output.VAL_HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        fastq_input = config['input']['fastq_pe'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        detection_method = config['detection_method'],
        citation_keys = config['citations']
    run:
        import datetime
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        # Add header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            params.sample_name,
            datetime.datetime.now(),
            params.pipeline_info['version'],
            ', '.join(entry['name'] for entry in params.fastq_input),
            [('Detection method', params.detection_method)],params.citation_keys['main']))

        # Add output sections
        report_structure = [
            ('Read trimming and basic QC', 'trim', [Path(input.report_downsampling), Path(input.report_trimming)]),
            ('Assembly', 'assem', [Path(input.report_assembly)]),
            ('Advanced QC', 'adv_qc', [Path(x) for x in (input.report_kraken, input.report_adv_qc)]),
            ('Variant calling', 'variant', [Path(input.report_variant)]),
            ('Resistance characterization', 'res', [Path(x) for x in (
                input.report_resfinder, input.report_argannot, input.report_card, input.report_ncbi_amr,
                input.report_pointfinder)]),
            ('Virulence characterization', 'viru', [Path(x) for x in (
                input.report_virulence, input.report_virulence_shiga)]),
            ('Subspecies identification', 'subspecies', [Path(x) for x in (
                input.report_species, input.report_subspecies, input.report_flexneri, input.report_subspecies_db)]),
            ('Plasmid replicon detection', 'plasmid', [Path(input.report_plasmidfinder)]),
            ('Sequence typing', 'st', [Path(x) for x in (
                input.report_mlst_warwick, input.report_mlst_pasteur, input.report_cgmlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ]
        SnakePipelineUtils.add_report_content(report, report_structure)


rule combine_summary_files:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        rules.init_summary.output.TSV,
        Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_SUMMARY,
        trimming.get_trimming_summary(config),
        Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_SUMMARY,
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_SUMMARY if 'kraken' in config['analyses'] else [],
        Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_SUMMARY,
        Path(config['working_dir']) / variant_filtering.OUTPUT_VARIANT_FILTERING_SUMMARY,
        # Gene detection
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='resfinder') if 'resfinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='card') if 'card' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='argannot') if 'argannot' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
        Path(config['working_dir']) / pointfinder.OUTPUT_POINTFINDER_SUMMARY if 'pointfinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='virulencefinder_shiga') if 'virulencefinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else [],
        # Sequence typing
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst_pasteur') if 'mlst_pasteur' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst_warwick') if 'mlst_warwick' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else []
    output:
        TSV = config.get('output_tabular')
    run:
        with open(output.TSV, 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
