from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import assembly_canu

camel = Camel.get_instance()


rule assembly_canu_run:
    """
    De-novo assembly using Canu.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        FASTA = Path(config['working_dir']) / 'assembly_canu' / 'canu' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'assembly_canu' / 'canu',
        canu_options = config.get('assembly', {}).get('canu', {}),
        read_type = 'SE'
    threads: 8
    priority: 1
    run:
        from camel.app.tools.canu.canu import Canu
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        canu = Canu(camel)

        # Reformat FASTQ dictionary
        fq_dict = SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_se='FASTQ', read_type=params.read_type)
        canu.add_input_files(fq_dict)
        step = Step(str(rule), canu, camel, params.running_dir)
        canu.update_parameters(**params.canu_options)
        canu.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(canu, output)

rule assembly_canu_filter_contig_length:
    """
    Filters out the small contigs.
    """
    input:
        FASTA = rules.assembly_canu_run.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / 'assembly_canu' / 'filtering' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_FILTERING_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'assembly_canu' / 'filtering',
        min_contig_length = config['assembly'].get('min_contig_length', 0) if 'assembly' in config else 0
    run:
        from camel.app.tools.seqtk.seqtkseq import SeqtkSeq
        seqtk = SeqtkSeq(camel)
        SnakemakeUtils.add_pickle_inputs(seqtk, input)
        step = Step(str(rule), seqtk, camel, params.running_dir)
        seqtk.update_parameters(output_filename='assembly_filtered.fasta', min_length=params.min_contig_length)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqtk, output)

rule assembly_canu_quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = rules.assembly_canu_filter_contig_length.output.FASTA
    output:
        TSV = Path(config['working_dir']) / 'assembly_canu' / 'quast' / 'tsv.io'
    params:
        running_dir = Path(config['working_dir']) / 'assembly_canu' / 'quast'
    run:
        from camel.app.tools.quast.quast import Quast
        quast = Quast(camel)
        SnakemakeUtils.add_pickle_inputs(quast, input)
        step = Step(str(rule), quast, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast, output)

rule assembly_canu_quast_extract_informs:
    """
    Extracts the information from the QUAST output file.
    """
    input:
        TSV = rules.assembly_canu_quast.output.TSV
    output:
        INFORMS = Path(config['working_dir']) / 'assembly_canu' / 'quast' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'assembly_canu' / 'quast'
    run:
        from camel.app.tools.quast.quastinformextractor import QuastInformExtractor
        quast_inform_extractor = QuastInformExtractor(camel)
        SnakemakeUtils.add_pickle_inputs(quast_inform_extractor, input)
        step = Step(str(rule), quast_inform_extractor, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast_inform_extractor, output)

rule assembly_canu_report:
    """
    Creates the HTML report for the assembly.
    """
    input:
        FASTA_Raw = rules.assembly_canu_run.output.FASTA,
        FASTA_Contig = rules.assembly_canu_filter_contig_length.output.FASTA,
        INFORMS_spades = rules.assembly_canu_run.output.INFORMS,
        INFORMS_quast = rules.assembly_canu_quast_extract_informs.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_REPORT
    params:
        running_dir = Path(config['working_dir']) / 'assembly_canu' / 'report',
        sample_name = config['sample_name']
    run:
        from camel.app.tools.pipelines.assembly.htmlreporterassembly import HtmlReporterAssembly
        from camel.app.io.tooliovalue import ToolIOValue
        reporter = HtmlReporterAssembly(camel)
        reporter.add_input_files({'SAMPLE_NAME': [ToolIOValue(params.sample_name)],
                                  'ASSEMBLER': [ToolIOValue('Canu')]})
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, camel, params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule assembly_canu_dump_summary_info:
    """
    Dumps the summary information from the assembly pipeline.
    """
    input:
        INFORMS_quast = rules.assembly_canu_quast_extract_informs.output.INFORMS
    output:
        Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_SUMMARY
    params:
        running_dir = Path(config['working_dir']) / 'assembly_canu' / 'summary'
    run:
        quast_informs = SnakemakeUtils.load_object(Path(input.INFORMS_quast))
        summary_data = [
            ('assembly_n50', quast_informs['contig']['N50']),
            ('assembly_nb_contigs', quast_informs['contig']['# contigs']),
            ('assembly_total_length', quast_informs['genome']['Total length'])
        ]
        with open(output[0], 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')

rule assembly_canu_map_reads:
    """
    Maps the reads to the canu assembly.
    """
    input:
        FQ = Path(config['working_dir']) / 'fq_dict.io',
        FASTA = rules.assembly_canu_filter_contig_length.output.FASTA
    output:
        SAM = Path(config['working_dir']) / 'assembly_canu' / 'minimap2' / 'sam.io'
    params:
        dir_ = Path(config['working_dir']) / 'assembly_canu' / 'minimap2'
    run:
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        minimap2 = Minimap2Mapping(Camel.get_instance())
        SnakemakeUtils.add_pickle_input(minimap2, 'FASTA', Path(input.FASTA))
        minimap2.add_input_files(SnakePipelineUtils.extracts_fq_input(
            Path(input.FQ), key_se='FASTQ', read_type='SE'))
        step = Step(str(rule), minimap2, Camel.get_instance(), params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(minimap2, output)

rule assembly_canu_sam_to_bam:
    """
    Converts the SAM file generated by bowtie2 to BAM format.
    """
    input:
        SAM = rules.assembly_canu_map_reads.output.SAM
    output:
        BAM = Path(config['working_dir']) / 'assembly_canu' / 'minimap2' / 'bam.io'
    params:
        running_dir = Path(config['working_dir']) / 'assembly_canu' / 'minimap2'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        samtools_view = SamtoolsView(camel)
        step = Step(str(rule), samtools_view, camel, params.running_dir)
        SnakemakeUtils.add_pickle_inputs(samtools_view, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_view, output)

rule assembly_canu_sort_bam:
    """
    Sorts the BAM alignment.
    """
    input:
        BAM = rules.assembly_canu_sam_to_bam.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'assembly_canu' / 'minimap2' / 'bam_sorted.io'
    params:
        running_dir = Path(config['working_dir']) / 'assembly_canu' / 'minimap2'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort(camel)
        step = Step(str(rule), samtools_sort, camel, params.running_dir)
        SnakemakeUtils.add_pickle_inputs(samtools_sort, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_sort, output)

rule assembly_canu_samtools_flagstat:
    """
    Runs samtools flagstat to determine the mapping rate.
    """
    input:
        BAM = rules.assembly_canu_sam_to_bam.output.BAM
    output:
        INFORMS = Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_MAPPING_RATE_INFORMS
    params:
        dir_ = Path(config['working_dir']) / 'assembly_canu' / 'samtools_flagstat'
    run:
        from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
        flagstat = SamtoolsFlagstat(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(flagstat, input)
        flagstat.run(Path(params.dir_).absolute())
        # Calculate mapping rate
        flagstat.informs['mapping_perc'] = 100 * flagstat.informs['mapped'][0] / flagstat.informs['total'][0]
        SnakemakeUtils.dump_tool_outputs(flagstat, output)

rule assembly_canu_samtools_depth:
    """
    Runs samtools depth on the BAM file of the reads mapped to the assembly.
    """
    input:
        BAM = rules.assembly_canu_sort_bam.output.BAM
    output:
        TSV = Path(config['working_dir']) / 'assembly_canu' / 'samtools_depth' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / assembly_canu.OUTPUT_ASSEMBLY_DEPTH_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'assembly_canu' / 'samtools_depth'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth(camel)
        step = Step(str(rule), samtools_depth, camel, params.running_dir)
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        step.run_step()
        samtools_depth.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)
