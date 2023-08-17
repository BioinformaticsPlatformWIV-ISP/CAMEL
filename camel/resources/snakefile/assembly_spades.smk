from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import assembly_spades

camel = Camel.get_instance()


rule assembly_spades_run:
    """
    De-novo assembly using SPAdes.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        FASTA_Contig = Path(config['working_dir']) / 'assembly_spades' / 'spades' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'assembly_spades' / 'spades',
        spades_options = config.get('assembly', {}).get('spades', {}),
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE'
    threads: 8
    priority: 1
    run:
        from camel.app.tools.spades.spades import SPAdes
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        spades = SPAdes(camel)

        # Reformat FASTQ dictionary
        fq_dict = SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE_1', keys_se=[
            'FASTQ_SE_1', 'FASTQ_SE_2'], key_se='FASTQ_SE_1', drop_empty=True, read_type=params.read_type)
        spades.add_input_files(fq_dict)
        step = Step(str(rule), spades, camel, params.running_dir)
        spades.update_parameters(**params.spades_options)
        spades.update_parameters(isolate=True, careful=False, threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spades, output)

rule assembly_filter_contig_length:
    """
    Filters out the small contigs.
    """
    input:
        FASTA = rules.assembly_spades_run.output.FASTA_Contig
    output:
        FASTA = Path(config['working_dir']) / 'assembly_spades' / 'filtering' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FILTERING_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'assembly_spades' / 'filtering',
        min_contig_length = config.get('assembly', {}).get('min_contig_length', 0)
    run:
        from camel.app.tools.seqtk.seqtkseq import SeqtkSeq
        seqtk = SeqtkSeq(camel)
        SnakemakeUtils.add_pickle_inputs(seqtk, input)
        step = Step(str(rule), seqtk, camel, params.running_dir)
        seqtk.update_parameters(output_filename='assembly_filtered.fasta', min_length=params.min_contig_length)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqtk, output)

rule assembly_quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = rules.assembly_filter_contig_length.output.FASTA
    output:
        TSV = Path(config['working_dir']) / 'assembly_spades' / 'quast' / 'tsv.io'
    params:
        running_dir = Path(config['working_dir']) / 'assembly_spades' / 'quast'
    run:
        from camel.app.tools.quast.quast import Quast
        quast = Quast(camel)
        SnakemakeUtils.add_pickle_inputs(quast, input)
        step = Step(str(rule), quast, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast, output)

rule assembly_quast_extract_informs:
    """
    Extracts the information from the QUAST output file.
    """
    input:
        TSV = rules.assembly_quast.output.TSV
    output:
        INFORMS = Path(config['working_dir']) / 'assembly_spades' / 'quast' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'assembly_spades' / 'quast'
    run:
        from camel.app.tools.quast.quastinformextractor import QuastInformExtractor
        quast_inform_extractor = QuastInformExtractor(camel)
        SnakemakeUtils.add_pickle_inputs(quast_inform_extractor, input)
        step = Step(str(rule), quast_inform_extractor, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast_inform_extractor, output)

rule assembly_report:
    """
    Creates the HTML report for the assembly.
    """
    input:
        FASTA_Raw = rules.assembly_spades_run.output.FASTA_Contig,
        FASTA_Contig = rules.assembly_filter_contig_length.output.FASTA,
        INFORMS_spades = rules.assembly_spades_run.output.INFORMS,
        INFORMS_quast = rules.assembly_quast_extract_informs.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'assembly_spades' / 'report',
        sample_name = config['sample_name']
    run:
        from camel.app.tools.pipelines.assembly.htmlreporterassembly import HtmlReporterAssembly
        from camel.app.io.tooliovalue import ToolIOValue
        reporter = HtmlReporterAssembly(camel)
        reporter.add_input_files({'SAMPLE_NAME': [ToolIOValue(params.sample_name)],
                                  'ASSEMBLER': [ToolIOValue('SPAdes')]})
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule assembly_dump_summary_info:
    """
    Dumps the summary information from the assembly pipeline.
    """
    input:
        INFORMS_spades = rules.assembly_spades_run.output.INFORMS,
        INFORMS_quast = rules.assembly_quast_extract_informs.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_SUMMARY
    params:
        running_dir = Path(config['working_dir']) / 'assembly_spades' / 'summary'
    run:
        spades_informs = SnakemakeUtils.load_object(Path(input.INFORMS_spades))
        quast_informs = SnakemakeUtils.load_object(Path(input.INFORMS_quast))
        summary_data = [
            ('assembly_n50', quast_informs['contig']['N50']),
            ('assembly_nb_contigs', quast_informs['contig']['# contigs']),
            ('assembly_total_length', quast_informs['genome']['Total length']),
            ('assembly_tool_version', spades_informs['_name'])
        ]
        with open(output.TSV, 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')

rule assembly_bt2_index:
    """
    Creates a bowtie2 index for the assembly.
    """
    input:
        FASTA_REF = rules.assembly_filter_contig_length.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'assembly_spades' / 'bowtie2' / 'genome_prefix.io'
    params:
        running_dir = Path(config['working_dir']) / 'assembly_spades' / 'bowtie2'
    run:
        from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index
        bowtie2_index = Bowtie2Index(camel)
        step = Step(str(rule), bowtie2_index, camel, params.running_dir)
        SnakemakeUtils.add_pickle_inputs(bowtie2_index, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bowtie2_index, output)

rule assembly_bt2_map:
    """
    Maps the reads against the assembled contigs.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
        INDEX_GENOME_PREFIX = rules.assembly_bt2_index.output.INDEX_GENOME_PREFIX
    output:
        BAM = Path(config['working_dir']) / 'assembly_spades' / 'bowtie2' / 'bam.io',
        INFORMS = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_MAPPING_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'assembly_spades' / 'bowtie2',
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE'
    threads: 8
    run:
        from camel.app.components.pipelines import pipeutils
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        # Bowtie 2
        bowtie2_map = Bowtie2Map(camel)
        bowtie2_map.add_input_files(SnakePipelineUtils.extracts_fq_input(
            Path(input.IO), key_se='FASTQ_SE', drop_empty=True, read_type=params.read_type))
        bowtie2_map.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(bowtie2_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))

        # Initialize tools
        samtools_view = SamtoolsView(Camel.get_instance())
        samtools_sort = SamtoolsSort(Camel.get_instance())
        samtools_sort.update_parameters(threads=threads)
        pipeutils.run_as_pipe([bowtie2_map, samtools_view, samtools_sort], Path(params.running_dir))

        # Save output
        SnakemakeUtils.dump_tool_output(samtools_sort, 'BAM', Path(output.BAM))
        bowtie2_map.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_object(bowtie2_map.informs, Path(output.INFORMS))

rule assembly_samtools_depth:
    """
    Runs samtools depth on the BAM file of the reads mapped to the assembly.
    """
    input:
        BAM = rules.assembly_bt2_map.output.BAM
    output:
        TSV = Path(config['working_dir']) / 'assembly_spades' / 'samtools_depth' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_DEPTH_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'assembly_spades' / 'samtools_depth'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth(camel)
        step = Step(str(rule), samtools_depth, camel, params.running_dir)
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        step.run_step()
        samtools_depth.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)
