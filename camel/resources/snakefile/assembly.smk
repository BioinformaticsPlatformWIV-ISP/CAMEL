from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import assembly, assembly_spades, assembly_flye


include: assembly_spades.SNAKEFILE_ASSEMBLY_SPADES
include: assembly_flye.SNAKEFILE_ASSEMBLY_FLYE


rule assembly_filter_contig_length:
    """
    Filters out the small contigs.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.get_fasta_raw(config)
    output:
        FASTA = Path(config['working_dir']) / 'assembly' / 'filtering' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / 'assembly' / 'filtering' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'assembly' / 'filtering',
        min_contig_length = config.get('assembly', {}).get('min_contig_length', 0)
    run:
        from camel.app.tools.seqtk.seqtkseq import SeqtkSeq
        seqtk = SeqtkSeq(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(seqtk, input)
        step = Step(str(rule), seqtk, Camel.get_instance(), params.dir_)
        seqtk.update_parameters(output_filename='assembly_filtered.fasta', min_length=params.min_contig_length)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqtk, output)

rule assembly_bt2_index:
    """
    Creates a bowtie2 index for the assembly.
    """
    input:
        FASTA_REF = rules.assembly_filter_contig_length.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'assembly' / 'bowtie2' / 'genome_prefix.io'
    params:
        dir_ = Path(config['working_dir']) / 'assembly' / 'bowtie2'
    run:
        from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index
        bowtie2_index = Bowtie2Index(Camel.get_instance())
        step = Step(str(rule), bowtie2_index, Camel.get_instance(), params.dir_)
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
        BAM = Path(config['working_dir']) / 'assembly' / 'bowtie2' / 'bam.io',
        INFORMS = Path(config['working_dir']) / 'assembly' / 'bowtie2' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'assembly' / 'bowtie2',
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE'
    threads: 8
    run:
        from camel.app.components.pipelines import pipeutils
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        # Bowtie 2
        bowtie2_map = Bowtie2Map(Camel.get_instance())
        bowtie2_map.add_input_files(SnakePipelineUtils.extracts_fq_input(
            Path(input.IO), key_se='FASTQ_SE', drop_empty=True, read_type=params.read_type))
        bowtie2_map.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(bowtie2_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))

        # Initialize tools
        samtools_view = SamtoolsView(Camel.get_instance())
        samtools_sort = SamtoolsSort(Camel.get_instance())
        samtools_sort.update_parameters(threads=threads)
        pipeutils.run_as_pipe([bowtie2_map, samtools_view, samtools_sort], Path(params.dir_))

        # Save output
        SnakemakeUtils.dump_tool_output(samtools_sort, 'BAM', Path(output.BAM))
        bowtie2_map.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_object(bowtie2_map.informs, Path(output.INFORMS))

rule assembly_minimap2_map:
    """
    Maps the reads to the Flye assembly using minimap2.
    """
    input:
        FQ = Path(config['working_dir']) / 'fq_dict.io',
        FASTA = rules.assembly_filter_contig_length.output.FASTA
    output:
        BAM = Path(config['working_dir']) / 'assembly' / 'minimap2' / 'bam.io',
        INFORMS = Path(config['working_dir']) / 'assembly' / 'minimap2' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'assembly' / 'minimap2'
    threads: 8
    run:
        from camel.app.components.pipelines import pipeutils
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        # Minimap2
        minimap2 = Minimap2Mapping(Camel.get_instance())
        SnakemakeUtils.add_pickle_input(minimap2, 'FASTA', Path(input.FASTA))
        minimap2.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.FQ), key_se='FASTQ', read_type='SE'))
        step = Step(str(rule), minimap2, Camel.get_instance(), params.dir_)
        step.run_step()

        # Initialize tools
        samtools_view = SamtoolsView(Camel.get_instance())
        samtools_sort = SamtoolsSort(Camel.get_instance())
        samtools_sort.update_parameters(threads=threads)
        pipeutils.run_as_pipe([minimap2, samtools_view, samtools_sort], Path(params.dir_))

        # Export output
        SnakemakeUtils.dump_tool_output(samtools_sort, 'BAM', Path(output.BAM))
        minimap2.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_object(minimap2.informs, Path(output.INFORMS))

rule assembly_samtools_depth:
    """
    Runs samtools depth on the BAM file of the reads mapped to the assembly.
    """
    input:
        BAM = Path(config['working_dir']) / 'assembly' / '{mapper}' / 'bam.io'
    output:
        TSV = Path(config['working_dir']) / 'assembly' / 'samtools_depth' / '{mapper}' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'assembly' / 'samtools_depth' / '{mapper}' / 'informs.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'assembly' / 'samtools_depth' / wildcards.mapper
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth(Camel.get_instance())
        step = Step(str(rule), samtools_depth, Camel.get_instance(), Path(str(params.dir_)))
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        step.run_step()
        samtools_depth.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)

rule assembly_samtools_flagstat:
    """
    Runs samtools flagstat to determine the mapping rate.
    """
    input:
        BAM = Path(config['working_dir']) / 'assembly' / '{mapper}' / 'bam.io'
    output:
        INFORMS = Path(config['working_dir']) / 'assembly' / 'samtools_flagstat' / '{mapper}' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'assembly' / 'samtools_flagstat'
    run:
        from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
        flagstat = SamtoolsFlagstat(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(flagstat, input)
        flagstat.run(Path(params.dir_).absolute())
        # Calculate mapping rate
        flagstat.informs['mapping_perc'] = 100 * flagstat.informs['mapped'][0] / flagstat.informs['total'][0]
        SnakemakeUtils.dump_tool_outputs(flagstat, output)

#################################################################
# The rules below generate a basic QUAST report / summary file. #
# For a more in-depth QUAST report use the quast.smk snakefile  #
#################################################################
rule assembly_quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = rules.assembly_filter_contig_length.output.FASTA
    output:
        TSV = Path(config['working_dir']) / 'assembly' / 'quast' / 'tsv.io'
    params:
        dir_ = Path(config['working_dir']) / 'assembly' / 'quast'
    run:
        from camel.app.tools.quast.quast import Quast
        quast = Quast(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(quast, input)
        step = Step(str(rule), quast, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast, output)

rule assembly_quast_extract_informs:
    """
    Extracts the information from the QUAST output file.
    """
    input:
        TSV = rules.assembly_quast.output.TSV
    output:
        INFORMS = Path(config['working_dir']) / 'assembly' / 'quast' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'assembly' / 'quast'
    run:
        from camel.app.tools.quast.quastinformextractor import QuastInformExtractor
        quast_inform_extractor = QuastInformExtractor(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(quast_inform_extractor, input)
        step = Step(str(rule), quast_inform_extractor, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast_inform_extractor, output)

rule assembly_quast_report:
    """
    Creates a basic HTML report for the assembly.
    """
    input:
        FASTA_raw = Path(config['working_dir']) / assembly.get_fasta_raw(config),
        FASTA_filt = rules.assembly_filter_contig_length.output.FASTA,
        INFORMS_assembler = assembly.get_command_informs(config),
        INFORMS_quast = rules.assembly_quast_extract_informs.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / 'assembly' / 'report' / 'html.io'
    params:
        dir_ = Path(config['working_dir']) / 'assembly' / 'report',
        sample_name = config['sample_name']
    run:
        from camel.app.tools.pipelines.assembly.htmlreporterassembly import HtmlReporterAssembly
        from camel.app.io.tooliovalue import ToolIOValue
        reporter = HtmlReporterAssembly(Camel.get_instance())
        reporter.add_input_files({'SAMPLE_NAME': [ToolIOValue(params.sample_name)]})
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule assembly_dump_summary_info:
    """
    Dumps the summary information from the assembly pipeline.
    """
    input:
        INFORMS_quast = rules.assembly_quast_extract_informs.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_SUMMARY
    run:
        quast_informs = SnakemakeUtils.load_object(Path(input.INFORMS_quast))
        summary_data = [
            ('assembly_n50', quast_informs['contig']['N50']),
            ('assembly_nb_contigs', quast_informs['contig']['# contigs']),
            ('assembly_total_length', quast_informs['genome']['Total length'])
        ]
        with open(output.TSV, 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')
