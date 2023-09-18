from pathlib import Path

from camel.app.camel import Camel
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import trimming, trimming_illumina, assembly_spades, quality_checks, \
    contamination_check_kraken, gene_detection, sequence_typing, downsampling, confindr, quast
from camel.scripts.neisseriapipeline.snakefile import serogroup_determination, gmats


#######################
# Included Snakefiles #
#######################
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: downsampling.SNAKEFILE_DOWNSAMPLING
include: confindr.SNAKEFILE_CONFINDR
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: assembly_spades.SNAKEFILE_ASSEMBLY_SPADES
include: quast.SNAKEFILE_QUAST
include: gene_detection.SNAKEFILE_GENE_DETECTION
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING
include: gmats.SNAKEFILE_GMATS
include: serogroup_determination.SNAKEFILE_SEROGROUP_DETERMINATION


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
        FASTA = Path(config['working_dir']) / gene_detection.INPUT_GENE_DETECTION_FASTA,
    shell:
        "cp {input.FASTA} {output.FASTA};"

rule link_fasta_to_typing:
    """
    This rule links the output of the assembly workflows to the sequence typing workflow.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA
    output:
        FASTA_typing = Path(config['working_dir']) / sequence_typing.INPUT_FASTA
    params:
        read_type = config['read_type']
    shell:
        """
        cp {input.FASTA} {output.FASTA_typing};
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

rule report_create_command_section:
    """
    Creates the report section containing the tool commands.
    """
    input:
        INFORMS_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_INFORMS,
        INFORMS_trimming = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_INFORMS,
        INFORMS_assembly = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_INFORMS,
        INFORMS_assembly_filt = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FILTERING_INFORMS,
        INFORMS_quast = Path(config['working_dir']) / quast.OUTPUT_QUAST_INFORMS,
        INFORMS_busco = Path(config['working_dir']) / quast.OUTPUT_BUSCO_INFORMS,
        INFORMS_kraken = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS if 'kraken' in config['analyses'] else [],
        INFORMS_confindr = Path(config['working_dir']) / confindr.OUTPUT_CONFINDR_INFORMS if 'confindr' in config['analyses'] else [],
        INFORMS_mapping = quality_checks.get_mapping_rate_informs(config),
        INFORMS_depth = quality_checks.get_depth_informs(config),
        INFORMS_resfinder = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='resfinder') if 'resfinder' in config['analyses'] else [],
        INFORMS_ncbi_amr = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
        INFORMS_mlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        INFORMS_rplf = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='rplf') if 'rplf' in config['analyses'] else [],
        INFORMS_bast = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='bast') if 'bast' in config['analyses'] else [],
        INFORMS_pora = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='pora') if 'pora' in config['analyses'] else [],
        INFORMS_porb = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='porb') if 'porb' in config['analyses'] else [],
        INFORMS_feta = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='feta') if 'feta' in config['analyses'] else [],
        INFORMS_amr = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='resistance_genes') if 'resistance_genes' in config['analyses'] else [],
        INFORMS_vaccine = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='vaccine_targets') if 'vaccine_targets' in config['analyses'] else [],
        INFORMS_fhbp = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='fhbp') if 'fhbp' in config['analyses'] else [],
        INFORMS_cgmlst = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        INFORMS_serogroup = Path(config['working_dir']) / serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_INFORMS if 'serogroup' in config['analyses'] else []
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        working_dir = config['working_dir']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        informs = []
        for content in [SnakemakeUtils.load_object(Path(io)) for io in input]:
            if type(content) is dict:
                informs.append(content)
            elif type(content) is list:
                informs.extend(content)
        section = SnakePipelineUtils.create_commands_section(informs, params.working_dir)
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule neisseria_additional_resistance_gene_metadata:
    """
    This rule is used to add resistance gene metadata for penA and rpoB genes.
    The data is parsed from the PubMLST webpage.
    """
    input:
        hits = Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_HITS).format(scheme='resistance_genes', locus_type='DNA', detection_method=config['detection_method']),
        VAL_HTML = sequence_typing.get_sequence_typing_report('resistance_genes', config),
        INFORMS_scheme = Path(config['working_dir']) / 'typing' / 'resistance_genes' / 'informs-locus_set.io'
    output:
        VAL_HTML = Path(config['working_dir']) / 'typing' / 'resistance_genes' / 'metadata' / 'report.html'
    params:
        working_dir = Path(config['working_dir']) / 'typing' / 'resistance_genes' / 'metadata',
        loci='penA, rpoB'
    run:
        from camel.app.pipeline.step import Step
        from camel.app.tools.pipelines.neisseria.resistancemetadataextractor import ResistanceMetadataExtractor
        extractor = ResistanceMetadataExtractor(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(extractor, input)
        step = Step(str(rule), extractor, Camel.get_instance(), params.working_dir, config)
        extractor.update_parameters(loci=params.loci)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(extractor, output)

rule combine_reports:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        report_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_REPORT,
        report_trimming = trimming.get_trimming_report(config),
        report_quast = Path(config['working_dir']) / quast.OUTPUT_QUAST_REPORT,
        report_kraken = Path(config['working_dir']) / (contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT if 'kraken' in config['analyses'] else contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY),
        report_confindr = Path(config['working_dir']) / (confindr.OUTPUT_CONFINDR_REPORT if 'confindr' in config['analyses'] else confindr.OUTPUT_CONFINDR_REPORT_EMPTY),
        report_adv_qc = Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_REPORT,
        report_resfinder = gene_detection.get_gene_detection_report('resfinder', config),
        report_ncbi_amr = gene_detection.get_gene_detection_report('ncbi_amr', config),
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        report_pora = sequence_typing.get_sequence_typing_report('pora', config),
        report_porb = sequence_typing.get_sequence_typing_report('porb', config),
        report_feta = sequence_typing.get_sequence_typing_report('feta', config),
        report_rplf = sequence_typing.get_sequence_typing_report('rplf', config),
        report_vaccine_targets = sequence_typing.get_sequence_typing_report('vaccine_targets', config),
        report_resistance_genes = rules.neisseria_additional_resistance_gene_metadata.output.VAL_HTML if 'resistance_genes' in config['analyses'] else Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_REPORT_EMPTY).format(scheme='resistance_genes'),
        report_fhbp = sequence_typing.get_sequence_typing_report('fhbp', config),
        report_bast = sequence_typing.get_sequence_typing_report('bast', config),
        report_gmats = Path(config['working_dir']) / (gmats.OUTPUT_GMATS_REPORT if 'gmats' in config['analyses'] else gmats.OUTPUT_GMATS_REPORT_EMPTY),
        report_serogroup = Path(config['working_dir']) / (serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_REPORT if 'serogroup' in config['analyses'] else serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_REPORT_EMPTY),
        report_serogroup_legacy = Path(config['working_dir']) / (serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_LEGACY_REPORT if 'serogroup' in config['analyses'] else serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_LEGACY_REPORT_EMPTY),
        report_citations = rules.report_pickle_citations.output.HTML,
        report_commands = rules.report_create_command_section.output.HTML
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
            params.pipeline_info['version'], ', '.join(entry['name'] for entry in params.fastq_input),
            [('Detection method', params.detection_method)], params.citation_keys['main']))

        # Add output sections
        report_structure = [
            ('Read trimming and basic QC', 'trim', [Path(input.report_downsampling), Path(input.report_trimming)]),
            ('Assembly', 'assem', [Path(input.report_quast)]),
            ('Advanced QC', 'adv_qc', [Path(x) for x in (
                input.report_kraken, input.report_confindr, input.report_adv_qc)]),
            ('Resistance characterization', 'res', [Path(x) for x in (input.report_resfinder, input.report_ncbi_amr)]),
            ('Sequence typing', 'st', [Path(x) for x in (
                input.report_mlst, input.report_rplf, input.report_pora, input.report_porb, input.report_feta,
                input.report_resistance_genes, input.report_vaccine_targets, input.report_fhbp, input.report_cgmlst)]),
            ('Antigen typing', 'at', [Path(x) for x in (input.report_bast, input.report_gmats)]),
            ('Serogroup determination', 'serogroup', [Path(
                input.report_serogroup), Path(input.report_serogroup_legacy)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ]
        SnakePipelineUtils.add_report_content(report, report_structure)

rule combine_summary_files:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        Path(config['working_dir']) / 'summary' / 'summary-init.tsv',
        Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_SUMMARY,
        trimming.get_trimming_summary(config),
        Path(config['working_dir']) / quast.OUTPUT_QUAST_SUMMARY,
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY,
        Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_SUMMARY if 'kraken' in config['analyses'] else [],
        Path(config['working_dir']) / confindr.OUTPUT_CONFINDR_SUMMARY if 'confindr' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='resfinder') if 'resfinder' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='rplf') if 'rplf' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='bast') if 'bast' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='pora') if 'pora' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='porb') if 'porb' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='feta') if 'feta' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='fhbp') if 'fhbp' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='resistance_genes') if 'resistance_genes' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='vaccine_targets') if 'vaccine_targets' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        Path(config['working_dir']) / gmats.OUTPUT_GMATS_SUMMARY if 'gmats' in config['analyses'] else [],
        Path(config['working_dir']) / serogroup_determination.OUTPUT_SEROGROUP_DETERMINATION_SUMMARY if 'serogroup' in config['analyses'] else []
    output:
        TSV = config.get('output_tabular')
    run:
        with open(output.TSV, 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
