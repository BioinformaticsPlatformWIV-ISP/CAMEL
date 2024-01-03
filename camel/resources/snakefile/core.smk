from pathlib import Path

from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import trimming_illumina, trimming_ont, assembly, gene_detection, downsampling, core, \
    sequence_typing, amrfinder, mobsuite

rule core_link_downsampling_input:
    """
    Creates the FASTQ input for the downsampling step. 
    """
    output:
        FASTQ_PE = Path(config['working_dir']) / str(downsampling.INPUT_DOWNSAMPLING_FASTQ).format(read_key='fastq_pe'),
        FASTQ_SE = Path(config['working_dir']) / str(downsampling.INPUT_DOWNSAMPLING_FASTQ).format(read_key='fastq_se')
    params:
        input_dict = config['input']
    run:
        from camel.app.io.tooliofile import ToolIOFile

        # FASTQ PE
        if 'fastq_pe' in params.input_dict:
            SnakemakeUtils.dump_object([ToolIOFile(Path(x['path'])) for x in params.input_dict['fastq_pe']], Path(output.FASTQ_PE))
        else:
            SnakemakeUtils.dump_object([], Path(output.FASTQ_PE))

        # FASTQ SE
        if 'fastq_se' in params.input_dict:
            SnakemakeUtils.dump_object([ToolIOFile(Path(params.input_dict['fastq_se'][0]['path']))], Path(output.FASTQ_SE))
        else:
            SnakemakeUtils.dump_object([], Path(output.FASTQ_SE))

rule core_link_trimmomatic_input:
    """
    Links the downsampling output to the input of the trimmomatic workflow.  
    """
    input:
        FASTQ = Path(config['working_dir']) / str(downsampling.OUTPUT_DOWNSAMPLING_FASTQ).format(read_key='fastq_pe')
    output:
        FASTQ = Path(config['working_dir']) / trimming_illumina.INPUT_TRIMMOMATIC_FASTQ
    shell:
        """
        cp {input.FASTQ} {output.FASTQ};
        """

rule core_link_ont_input:
    """
    Links the downsampling output to the input of the trimmomatic workflow.  
    """
    input:
        FASTQ = Path(config['working_dir']) / str(downsampling.OUTPUT_DOWNSAMPLING_FASTQ).format(read_key='fastq_se')
    output:
        FASTQ = Path(config['working_dir']) / trimming_ont.INPUT_ONT_FASTQ
    shell:
        """
        cp {input.FASTQ} {output.FASTQ};
        """

rule core_select_fasta:
    """
    This rules links the output of the assembly workflow to the other workflows. 
    """
    input:
        FASTA = Path(config['working_dir'], assembly.OUTPUT_ASSEMBLY_FASTA)
    output:
        FASTA = Path(config['working_dir']) / gene_detection.INPUT_GENE_DETECTION_FASTA
    shell:
        "cp {input.FASTA} {output.FASTA};"

rule core_link_fasta_to_typing:
    """
    This rule core_links the output of the assembly workflows to the sequence typing workflow.
    """
    input:
        FASTA = Path(config['working_dir'], assembly.OUTPUT_ASSEMBLY_FASTA)
    output:
        FASTA_typing = Path(config['working_dir']) / sequence_typing.INPUT_FASTA
    shell:
        """
        cp {input.FASTA} {output.FASTA_typing};
        """

rule core_link_fasta_to_amrfinder:
    """
    This rule core_links the output of the assembly workflows to the AMRFinder workflow.
    """
    input:
        FASTA = Path(config['working_dir'], assembly.OUTPUT_ASSEMBLY_FASTA)
    output:
        FASTA = Path(config['working_dir']) / amrfinder.INPUT_AMRFINDER_FASTA
    shell:
        """
        cp {input.FASTA} {output.FASTA};
        """

rule core_link_fasta_to_mob_suite:
    """
    This rule core_links the output of the assembly workflows to the MOB-suite workflow.
    """
    input:
        FASTA = Path(config['working_dir'], assembly.OUTPUT_ASSEMBLY_FASTA)
    output:
        FASTA = Path(config['working_dir']) / mobsuite.INPUT_MOBSUITE_FASTA
    shell:
        """
        cp {input.FASTA} {output.FASTA};
        """

rule core_collect_trimmed_fastq_data:
    """
    This rule core_creates an IO object with the trimmed FASTQ files.
    Other workflows such as 'Kraken' or 'Assembly' rely on this dictionary to get input files (PE or SE).
    """
    input:
        unpack(lambda _: core.get_fq_input(config['input_type'], config['working_dir']))
    output:
        IO_FASTQ = Path(config['working_dir']) / 'fq_dict.io'
    params:
        input_type = config['input_type']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.construct_fq_dict(input, params.input_type, Path(output.IO_FASTQ))

rule core_init_summary:
    """
    Initializes the summary output file.
    """
    output:
        TSV = Path(config['working_dir']) / 'summary' / 'summary-init.tsv'
    run:
        import datetime
        from camel.app.components.pipelines.reportpipeline import ReportPipeline

        analysis_date = datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT)
        content = [
            ('pipeline_name', config['pipeline']['name']),
            ('pipeline_version', config['pipeline']['version']),
            ('input_type', config['input_type']),
            ('sample', config['sample_name']),
            ('input_files', ReportPipeline.format_input_string(config['input'])),
            ('analysis_date', analysis_date),
            ('detection_method', config['detection_method'])
        ]
        with open(output.TSV, 'w') as handle:
            for kv_pair in content:
                handle.write('\t'.join(kv_pair))
                handle.write('\n')

rule core_report_pickle_citations:
    """
    This rule core_creates a pickle with a report section containing the citations.
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
