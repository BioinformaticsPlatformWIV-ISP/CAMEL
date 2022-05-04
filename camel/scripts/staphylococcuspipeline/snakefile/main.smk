from pathlib import Path

from camel.resources.snakefile import trimming_illumina, assembly_spades, gene_detection, trimming, \
    contamination_check_kraken, quality_checks, sequence_typing, pointfinder, plasmidspades, lrefinder, downsampling
from camel.scripts.staphylococcuspipeline.snakefile import spatyping, sccmectyping

#######################
# Included Snakefiles #
#######################
include: downsampling.SNAKEFILE_DOWNSAMPLING
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: assembly_spades.SNAKEFILE_ASSEMBLY_SPADES
include: pointfinder.SNAKEFILE_POINTFINDER
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING
include: spatyping.SNAKEFILE_SPATYPING
include: sccmectyping.SNAKEFILE_SCCMEC_TYPING
include: plasmidspades.SNAKEFILE_PLASMID_SPADES
include: lrefinder.SNAKEFILE_LREFINDER


rule all:
    """
    Rule to generate the required output files.
    """
    input:
        HTML = config['output_report'],
        TSV = config['output_tabular']

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
    Links the downsmapling output to the input of the trimmomatic workflow.  
    """
    input:
        FASTQ = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_FASTQ
    output:
        FASTQ = Path(config['working_dir']) / trimming_illumina.INPUT_TRIMMOMATIC_FASTQ
    shell:
        """
        cp {input.FASTQ} {output.FASTQ};
        """

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
    """
    Creates a report section with the commands used in the pipeline. 
    """
    input:
        INFORMS_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_INFORMS,
        INFORMS_trimming = trimming.get_trimming_command_informs(config),
        INFORMS_assembly = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_INFORMS,
        INFORMS_assembly_filt = Path(config['working_dir']) / 'assembly_spades' / 'filtering' / 'informs.io',
        INFORMS_kraken = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS if 'kraken' in config['analyses'] else [],
        INFORMS_mapping = quality_checks.get_mapping_rate_informs(config),
        INFORMS_depth = quality_checks.get_depth_informs(config),
        INFORMS_resfinder = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='resfinder') if 'resfinder' in config['analyses'] else [],
        INFORMS_ncbi_amr = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
        INFORMS_lrefinder =Path(config['working_dir']) / lrefinder.OUTPUT_LREFINDER_INFORMS if 'lrefinder' in config['analyses'] else [],
        INFORMS_vfdb_core = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        INFORMS_se_toxins = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='se_toxins') if 'se_toxins' in config['analyses'] else [],
        INFORMS_pointfinder = Path(config['working_dir']) / pointfinder.OUTPUT_POINTFINDER_INFORMS if 'pointfinder' in config['analyses'] else [],
        INFORMS_plasmidspades =Path(config['working_dir']) / plasmidspades.OUTPUT_PLASMIDSPADES_INFORMS if 'plasmidspades' in config['analyses'] else []
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
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

rule report_combine_all:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        report_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_REPORT,
        report_trimming = trimming.get_trimming_report(config),
        report_assembly = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_REPORT,
        report_kraken = Path(config['working_dir']) / (contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT if 'kraken' in config['analyses'] else contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY),
        report_adv_qc = Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT,
        # spa typing
        report_spa_typing = Path(config['working_dir']) / (spatyping.OUTPUT_SPATYPING_REPORT if 'spa_typing' in config['analyses'] else spatyping.OUTPUT_SPATYPING_REPORT_EMPTY),
        # SCCmec typing
        report_sccmec_genes = gene_detection.get_gene_detection_report('sccmec_genes', config, 'sccmec_typing'),
        report_sccmec_cassette = gene_detection.get_gene_detection_report('sccmec_cassette', config, 'sccmec_typing'),
        report_sccmec_typing = Path(config['working_dir']) / (sccmectyping.OUTPUT_SCCMEC_TYPING_REPORT if 'sccmec_typing' in config['analyses'] else sccmectyping.OUTPUT_SCCMEC_TYPING_REPORT_EMPTY),
        # AMR detection
        report_resfinder = gene_detection.get_gene_detection_report('resfinder', config),
        report_ncbi_amr = gene_detection.get_gene_detection_report('ncbi_amr', config),
        report_pointfinder = Path(config['working_dir']) / (pointfinder.OUTPUT_POINTFINDER_REPORT if 'pointfinder' in config['analyses'] else pointfinder.OUTPUT_POINTFINDER_REPORT_EMPTY),
        report_lrefinder = Path(config['working_dir']) / (lrefinder.OUTPUT_LREFINDER_REPORT if 'lrefinder' in config['analyses'] else lrefinder.OUTPUT_LREFINDER_REPORT_EMPTY),
        # Virulence detection
        report_vf_exoenzyme = gene_detection.get_gene_detection_report('vf_exoenzyme', config, 'virulencefinder'),
        report_vf_hostimm = gene_detection.get_gene_detection_report('vf_hostimm', config, 'virulencefinder'),
        report_vf_toxin = gene_detection.get_gene_detection_report('vf_toxin', config, 'virulencefinder'),
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        report_se_toxins = gene_detection.get_gene_detection_report('se_toxins', config),
        # Plasmid characterization
        report_plasmidfinder = gene_detection.get_gene_detection_report('plasmidfinder', config),
        report_plasmidspades = Path(config['working_dir']) / (
            plasmidspades.OUTPUT_PLASMIDSPADES_REPORT if 'plasmidspades' in config['analyses'] else plasmidspades.OUTPUT_PLASMIDSPADES_REPORT_EMPTY),
        report_plasmid_resfinder = Path(config['working_dir']) / str(
            plasmidspades.OUTPUT_PLASMIDSPADES_GENE_DETECTION_REPORT if 'plasmidspades' in config['analyses'] else plasmidspades.OUTPUT_PLASMIDSPADES_GENE_DETECTION_REPORT_EMPTY).format(db='resfinder'),
        report_plasmid_ndaro = Path(config['working_dir']) / str(
            plasmidspades.OUTPUT_PLASMIDSPADES_GENE_DETECTION_REPORT if 'plasmidspades' in config['analyses'] else plasmidspades.OUTPUT_PLASMIDSPADES_GENE_DETECTION_REPORT_EMPTY).format(db='ncbi_amr'),
        # Typing
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        # Report
        report_citations = rules.report_pickle_citations.output.HTML,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        fastq_input = config['input']['fastq_pe'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        detection_method = config['detection_method']
    run:
        import datetime
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        # Add header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            params.sample_name, datetime.datetime.now(), params.pipeline_info['version'],
            ', '.join(entry['name'] for entry in params.fastq_input), [('Detection method', params.detection_method)]))

        # Add report content
        report_structure = [
            ('Read trimming and basic QC', 'trim', [Path(input.report_downsampling), Path(input.report_trimming)]),
            ('Assembly', 'assem', [Path(input.report_assembly)]),
            ('Advanced QC', 'adv_qc', [Path(x) for x in (input.report_kraken, input.report_adv_qc)]),
            ('<i>spa</i> typing', 'spa', [Path(input.report_spa_typing)]),
            ('SCC<i>mec</i> typing', 'sccmec', [Path(x) for x in (
                input.report_sccmec_genes, input.report_sccmec_cassette, input.report_sccmec_typing)]),
            ('AMR detection', 'amr', [Path(x) for x in (
                input.report_resfinder, input.report_ncbi_amr, input.report_pointfinder, input.report_lrefinder)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (
                input.report_vfdb_core, input.report_vf_exoenzyme, input.report_vf_hostimm, input.report_vf_toxin,
                input.report_se_toxins)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_plasmidfinder, input.report_plasmidspades, input.report_plasmid_resfinder,
                input.report_plasmid_ndaro)]),
            ('Sequence typing', 'st', [Path(x) for x in (input.report_mlst, input.report_cgmlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ]
        SnakePipelineUtils.add_report_content(report, report_structure)

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
    Combines the summary information of several steps into a single TSV file.
    """
    input:
        rules.summary_init.output.TSV,
        Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_SUMMARY,
        trimming.get_trimming_summary(config),
        Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_SUMMARY,
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_SUMMARY if 'kraken' in config['analyses'] else [],
        # spa and SCCmec typing
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='sccmec_genes') if 'sccmec_typing' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='sccmec_cassette') if 'sccmec_typing' in config['analyses'] else [],
        Path(config['working_dir']) / spatyping.OUTPUT_SPATYPING_SUMMARY if 'spa_typing' in config['analyses'] else [],
        Path(config['working_dir']) / sccmectyping.OUTPUT_SCCMEC_TYPING_SUMMARY if 'sccmec_typing' in config['analyses'] else [],
        # AMR detection
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='resfinder') if 'resfinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
        Path(config['working_dir']) / pointfinder.OUTPUT_POINTFINDER_SUMMARY if 'pointfinder' in config['analyses'] else [],
        Path(config['working_dir']) / lrefinder.OUTPUT_LREFINDER_SUMMARY if 'lrefinder' in config['analyses'] else [],
        # Virulence detection
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='vf_exoenzyme') if 'virulencefinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='vf_hostimm') if 'virulencefinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='vf_toxin') if 'virulencefinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='se_toxins') if 'se_toxins' in config['analyses'] else [],
        # Plasmid characterization
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else [],
        Path(config['working_dir']) / plasmidspades.OUTPUT_PLASMIDSPADES_SUMMARY if 'plasmidspades' in config['analyses'] else [],
        # Sequence typing
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else []
    output:
        TSV = config['output_tabular']
    run:
        with open(output[0], 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
