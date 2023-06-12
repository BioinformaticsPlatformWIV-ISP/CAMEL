import shutil
from pathlib import Path

import pandas as pd

from camel.resources.snakefile import trimming, trimming_illumina, assembly_spades, assembly_canu, \
    quality_checks, contamination_check_kraken, sequence_typing, downsampling, amrfinder, trimming_ont, gene_detection, \
    mobsuite
from camel.scripts.bacilluspipeline.snakefile import btyper, ani

#######################
# Included Snakefiles #
#######################
include: downsampling.SNAKEFILE_DOWNSAMPLING
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: trimming_ont.SNAKEFILE_TRIMMING_ONT
include: contamination_check_kraken.SNAKEFILE_CONTAMINATION_CHECK_KRAKEN
include: quality_checks.SNAKEFILE_QUALITY_CHECKS
include: assembly_spades.SNAKEFILE_ASSEMBLY_SPADES
include: assembly_canu.SNAKEFILE_ASSEMBLY_CANU
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

rule link_downsampling_input:
    """
    Creates the FASTQ input for the downsampling step.
    """
    output:
        FASTQ = Path(config['working_dir']) / downsampling.INPUT_DOWNSAMPLING_FASTQ
    params:
        config_input=config['input'],
        read_type=config['read_type']
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        if params.read_type == 'illumina':
            SnakemakeUtils.dump_object([ToolIOFile(Path(x['path'])) for x in config['input']['fastq_pe']],
                Path(output.FASTQ))
        elif params.read_type == 'nanopore':
            SnakemakeUtils.dump_object([ToolIOFile(Path(x['path'])) for x in config['input']['fastq_se']],
                Path(output.FASTQ))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

rule link_downsampling_to_trimming_workflows:
    """
    Links the downsampling output to the input of the ONT trimming workflow.  
    """
    input:
        FASTQ = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_FASTQ
    output:
        FASTQ_ilmn = Path(config['working_dir']) / trimming_illumina.INPUT_TRIMMOMATIC_FASTQ if config['read_type'] == 'illumina' else [],
        FASTQ_ont = Path(config['working_dir']) / trimming_ont.INPUT_ONT_FASTQ if config['read_type'] == 'nanopore' else []
    params:
        read_type=config['read_type']
    run:
        if params.read_type == 'nanopore':
            shutil.copyfile(Path(input.FASTQ),Path(output.FASTQ_ont))
        elif params.read_type == 'illumina':
            shutil.copyfile(Path(input.FASTQ),Path(output.FASTQ_ilmn))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

rule select_fastq_to_io:
    """
    Creates an IO object with the trimmed FASTQ files.
    Other workflows such as Kraken or de novo assembly rely on this dictionary to get input files (PE or SE).
    """
    input:
        FASTQ_ilmn = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_DICT if config['read_type'] == 'illumina' else [],
        FASTQ_ont = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_DICT if config['read_type'] == 'nanopore' else []
    output:
        IO_FASTQ = Path(config['working_dir']) / 'fq_dict.io'
    params:
        read_type = config['read_type']
    run:
        if params.read_type == 'illumina':
            shutil.copyfile(Path(input.FASTQ_ilmn), Path(output.IO_FASTQ))
        elif params.read_type == 'nanopore':
            shutil.copyfile(Path(input.FASTQ_ont), Path(output.IO_FASTQ))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

rule select_fasta_to_gene_detection:
    """
    This rule links the output of the assembly workflows to the gene detection workflow.
    """
    input:
        FASTA_spades = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA if config['read_type'] == 'illumina' else [],
        FASTA_canu = Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_FASTA if config['read_type'] == 'nanopore' else []
    output:
        FASTA_genedetection = Path(config['working_dir']) / gene_detection.INPUT_GENE_DETECTION_FASTA
    params:
        read_type = config['read_type']
    run:
        if params.read_type == 'nanopore':
            shutil.copyfile(Path(input.FASTA_canu), Path(output.FASTA_genedetection))
        elif params.read_type == 'illumina':
            shutil.copyfile(Path(input.FASTA_spades), Path(output.FASTA_genedetection))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

rule select_fasta_to_tools:
    """
    This rule links the output of the assembly workflow to the amrfinder and mobsuite workflows.
    """
    input:
        FASTA_spades = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA if config['read_type'] == 'illumina' else [],
        FASTA_canu = Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_FASTA if config['read_type'] == 'nanopore' else []
    output:
        FASTA_amrfinder = Path(config['working_dir']) / amrfinder.INPUT_AMRFINDER_FASTA,
        FASTA_mobsuite = Path(config['working_dir']) / mobsuite.INPUT_MOBSUITE_FASTA
    params:
        read_type=config['read_type']
    run:
        if params.read_type == 'nanopore':
            shutil.copyfile(Path(input.FASTA_canu), Path(output.FASTA_amrfinder))
            shutil.copyfile(Path(input.FASTA_canu), Path(output.FASTA_mobsuite))
        elif params.read_type == 'illumina':
            shutil.copyfile(Path(input.FASTA_spades), Path(output.FASTA_amrfinder))
            shutil.copyfile(Path(input.FASTA_spades), Path(output.FASTA_mobsuite))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

rule select_fasta_to_tools_subtilis:
    """
    This rule links the output of the assembly workflow to the fastANI workflow if the species is B. subtilis.
    """
    input:
        FASTA_spades = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA if config['read_type'] == 'illumina' else [],
        FASTA_canu = Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_FASTA if config['read_type'] == 'nanopore' else []
    output:
        FASTA_ani = Path(config['working_dir']) / ani.INPUT_FASTA_ANI
    params:
        read_type = config['read_type']
    run:
        if params.read_type == 'nanopore':
            shutil.copyfile(Path(input.FASTA_canu), Path(output.FASTA_ani))
        elif params.read_type == 'illumina':
            shutil.copyfile(Path(input.FASTA_spades), Path(output.FASTA_ani))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

rule select_fasta_to_tools_cereus:
    """
    This rule links the output of the assembly workflow to the BTyper workflow if the species is B. cereus.
    """
    input:
        FASTA_spades = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA if config['read_type'] == 'illumina' else [],
        FASTA_canu = Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_FASTA if config['read_type'] == 'nanopore' else []
    output:
        FASTA_btyper = Path(config['working_dir']) / btyper.INPUT_BTYPER_FASTA
    params:
        read_type=config['read_type']
    run:
        if params.read_type == 'nanopore':
            shutil.copyfile(Path(input.FASTA_canu), Path(output.FASTA_btyper))
        elif params.read_type == 'illumina':
            shutil.copyfile(Path(input.FASTA_spades),Path(output.FASTA_btyper))
        else:
            raise ValueError(f'Unsupported read type: {params.read_type}')

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
    This rule retrieves the commands for the final report.
    """
    input:
        INFORMS_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_INFORMS,
        INFORMS_trimming = trimming.get_trimming_command_informs(config),
        INFORMS_assembly = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_INFORMS if config['read_type'] == 'illumina' else Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_INFORMS,
        INFORMS_assembly_filt = Path(config['working_dir']) / 'assembly_spades' / 'filtering' / 'informs.io' if config['read_type'] == 'illumina' else [],
        INFORMS_kraken = Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS if 'kraken' in config['analyses'] else [],
        INFORMS_gmo = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='gmo') if 'gmo' in config['analyses'] else [],
        INFORMS_vfdb_core = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        INFORMS_mapping = quality_checks.get_mapping_rate_informs(config) if config['read_type'] == 'illumina' else [],
        INFORMS_depth = quality_checks.get_depth_informs(config) if config['read_type'] == 'illumina' else [],
        INFORMS_btyper = Path(config['working_dir']) / btyper.OUTPUT_INFORMS_BTYPER if config['contamination_check']['expected_species'] == 'Bacillus cereus' else [],
        INFORMS_ani = Path(config['working_dir']) / str(ani.OUTPUT_INFORMS_ANI) if config['contamination_check']['expected_species'] == 'Bacillus subtilis' else [],
        INFORMS_amrfinder = Path(config['working_dir']) / str(amrfinder.OUTPUT_AMRFINDER_INFORMS) if 'amrfinder' in config['analyses'] else [],
        INFORMS_mob_suite = Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_INFORMS if 'mob_suite' in config['analyses'] else []
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
            [('Detection method', params.detection_method)],
            params.citation_keys['main']
        )
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.HTML))

checkpoint determine_species:
    input:
        TSV = Path(config['working_dir']) / 'contamination_check' / 'kraken2' / 'tsv-report.io'
    output:
        TXT = Path(config['working_dir']) / 'detected_species.txt'
    run:
        path_tsv = SnakemakeUtils.load_object(Path(input.TSV))[0].path
        data_k2 = pd.read_table(path_tsv, usecols=[0, 3, 5], names=['perc', 'level', 'name'])
        species = data_k2[data_k2['level'] == 'S'].sort_values(by='perc', ascending=False).iloc[0]['name'].strip()
        with open(output.TXT, 'w') as handle:
            handle.write(species)
            handle.write('\n')

rule report_content_cereus:
    """
    Creates the main content of the report when the detected species is Bacillus cereus. 
    """
    input:
        report_init = rules.report_init.output.HTML,
        report_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_REPORT,
        report_trimming = trimming.get_trimming_report(config),
        report_btyper=Path(config['working_dir']) / (btyper.OUTPUT_BTYPER_REPORT if 'btyper' in config['analyses'] else btyper.OUTPUT_BTYPER_REPORT_EMPTY),
        report_vfdb_core=gene_detection.get_gene_detection_report('vfdb_core', config)
    output:
        HTML = Path(config['working_dir']) / 'report' / 'report_cereus.html'
    params:
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline']
    run:
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakemakeUtils.load_object(Path(input.report_init))[0].value)
        report_structure = [
            ('Read trimming and basic QC', 'trim', [Path(input.report_downsampling), Path(input.report_trimming)]),
            ('BTyper results', 'btyper', [Path(input.report_btyper)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (input.report_vfdb_core,)])
        ]
        SnakePipelineUtils.add_report_content(report, report_structure)
        report.save()

rule report_content_subtilis:
    """
    Creates the main content of the report when the detected species is Bacillus subtilis. 
    """
    input:
        report_init = rules.report_init.output.HTML,
        report_downsampling = Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_REPORT,
        report_trimming = trimming.get_trimming_report(config),
        report_amrfinder= Path(config['working_dir']) / amrfinder.OUTPUT_AMRFINDER_REPORT,
        report_vfdb_core= gene_detection.get_gene_detection_report('vfdb_core', config)

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
            ('Read trimming and basic QC', 'trim', [Path(input.report_downsampling), Path(input.report_trimming)]),
            ('Antimicrobial resistance detection', 'amr', [Path(input.report_amrfinder)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (input.report_vfdb_core,)])
        ]
        SnakePipelineUtils.add_report_content(report,report_structure)
        report.save()

def select_report_input(wildcards) -> Path:
    with open(checkpoints.determine_species.get().output.TXT) as h:
        species_k2 = h.readline().strip()
    if species_k2 == 'Bacillus subtilis':
        return Path(config['working_dir']) / 'report' / 'report_subtilis.html'
    elif species_k2 == 'Bacillus cereus':
        return Path(config['working_dir']) / 'report' / 'report_cereus.html'
    else:
        raise ValueError(f'Unsupported species: {species_k2}')

rule report_select_by_species:
    """
    Selects the report content based on the detected species.
    """
    input:
        HTML = select_report_input
    output:
        HTML = config['output_report']
    params:
        output_dir = config['output_dir']
    shell:
        """
        cp {input.HTML} {output.HTML};
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
        Path(config['working_dir']) / downsampling.OUTPUT_DOWNSAMPLING_SUMMARY,
        trimming.get_trimming_summary(config),
        Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_SUMMARY if config['read_type'] == 'illumina' else Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_SUMMARY,
        Path(config['working_dir']) / quality_checks.OUTPUT_QUALITY_CHECKS_SUMMARY if config['read_type'] == 'illumina' else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='gmo') if 'gmo' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_SUMMARY).format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else [],
        Path(config['working_dir']) / contamination_check_kraken.OUTPUT_CONTAMINATION_SUMMARY if 'kraken' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_SUMMARY).format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        Path(config['working_dir']) / btyper.OUTPUT_BTYPER_SUMMARY if config['contamination_check']['expected_species'] == 'Bacillus cereus' else [],
        Path(config['working_dir']) / ani.OUTPUT_ANI_SUMMARY if config['contamination_check']['expected_species'] == 'Bacillus subtilis' else [],
        Path(config['working_dir']) / amrfinder.OUTPUT_AMRFINDER_SUMMARY if 'amrfinder' in config['analyses'] else [],
        Path(config['working_dir']) / mobsuite.OUTPUT_MOB_SUITE_SUMMARY if 'mob_suite' in config['analyses'] else []
    output:
        config.get('output_tabular')
    run:
        with open(output[0], 'w') as handle_out:
            for summary_input in input:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
