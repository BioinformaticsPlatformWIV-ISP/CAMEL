import shutil
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import short_read_polishing

camel = Camel.get_instance()


rule polishing_copy_fasta:
    """
    Copies the input FASTA file into the polypolish input folder.
    """
    input:
        FASTA = lambda wildcards: Path(config['working_dir']) / str(short_read_polishing.INPUT_ASSEMBLY_FASTA).format(assembly_type=wildcards.assembly_type)
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'polypolish' / 'input_genome.fasta'
    run:
        fasta = SnakemakeUtils.load_object(Path(str(input.FASTA)))
        fasta_file = fasta[0].path
        shutil.copyfile(fasta_file, output.FASTA)

rule polishing_samtools_index_polypolish:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA_REF = rules.polishing_copy_fasta.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'polypolish' / 'input_genome.fasta.fai'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'polypolish'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex(camel)
        samtools.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule), samtools, camel, Path(str(params.running_dir)))
        step.run_step()

rule polishing_bwa_index:
    """
    Creates a bwa index for the assembly.
    """
    input:
        FASTA_REF = rules.polishing_copy_fasta.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'polypolish' / 'genome_prefix.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'polypolish'
    run:
        from camel.app.tools.bwa.bwaindex import BWAIndex
        bwa = BWAIndex(camel)
        bwa.add_input_files({'FASTA_REF': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule), bwa, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bwa, output)

rule polishing_read_mapping_1:
    """
    Maps the forward reads against the assembly.
    """
    input:
        FQ_dict = Path(config['working_dir']) / 'fq_dict.io',
        INDEX_GENOME_PREFIX_BWA = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX,
        INDEX_GENOME_PREFIX_SAMTOOLS = rules.polishing_samtools_index_polypolish.output.INDEX_GENOME_PREFIX,
        FASTA = rules.polishing_copy_fasta.output.FASTA,
        INDEX_GENOME_PREFIX = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'read_mapping' / 'forward' / 'bwa_readmap.sam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'read_mapping' / 'forward'
    threads: 8
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.bwa.bwamap import BWAMap
        bwa_map = BWAMap(camel)
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict), 'illumina')
        bwa_map.add_input_files({'FASTQ_SE': [fq_in.pe[0]]})
        bwa_map.update_parameters(threads=threads, all_alns=True)
        SnakemakeUtils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step = Step(str(rule), bwa_map, camel, Path(str(params.running_dir)))
        step.run_step()

rule polishing_read_mapping_2:
    """
    Maps the forward reads against the assembly.
    """
    input:
        FQ_dict = Path(config['working_dir']) / 'fq_dict.io',
        INDEX_GENOME_PREFIX_BWA = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX,
        INDEX_GENOME_PREFIX_SAMTOOLS = rules.polishing_samtools_index_polypolish.output.INDEX_GENOME_PREFIX,
        FASTA = rules.polishing_copy_fasta.output.FASTA,
        INDEX_GENOME_PREFIX = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'read_mapping' / 'reverse' / 'bwa_readmap.sam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'read_mapping' / 'reverse'
    threads: 8
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.bwa.bwamap import BWAMap
        bwa_map = BWAMap(camel)
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict), 'illumina')
        bwa_map.add_input_files({'FASTQ_SE': [fq_in.pe[1]]})
        bwa_map.update_parameters(threads=threads, all_alns=True)
        SnakemakeUtils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step = Step(str(rule), bwa_map, camel, Path(str(params.running_dir)))
        step.run_step()

rule polishing_polypolish_insert_filter:
    input:
        SAM_1 = rules.polishing_read_mapping_1.output.SAM,
        SAM_2 = rules.polishing_read_mapping_2.output.SAM
    output:
        SAM_1 = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'read_mapping' / 'alignment_filtered_1.sam',
        SAM_2 = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'read_mapping' / 'alignment_filtered_2.sam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'read_mapping'
    threads: 8
    run:
        from camel.app.tools.polypolish.polypolishinsertfilter import PolypolishInsertFilter
        insert_filter = PolypolishInsertFilter(camel)
        insert_filter.add_input_files({'SAM': [ToolIOFile(Path(input.SAM_1)), ToolIOFile(Path(input.SAM_2))]})
        step = Step(str(rule), insert_filter, camel, Path(str(params.running_dir)))
        step.run_step()

rule polishing_polypolish:
    """
    First polishing with polypolish.
    """
    input:
        SAM_1 = rules.polishing_polypolish_insert_filter.output.SAM_1,
        SAM_2 = rules.polishing_polypolish_insert_filter.output.SAM_2,
        FASTA = rules.polishing_copy_fasta.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'polypolish' / 'polished.fasta',
        INFORMS = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'polypolish'  / 'polypolish.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'polypolish',
        polypolish_options = config.get('polishing', {}).get('polypolish', {})
    run:
        from camel.app.tools.polypolish.polypolish import Polypolish
        polypolish = Polypolish(camel)
        polypolish.add_input_files({
            'SAM': [ToolIOFile(Path(input.SAM_1)), ToolIOFile(Path(input.SAM_2))],
            'FASTA':[ToolIOFile(Path(input.FASTA))]})
        polypolish.update_parameters(**params.polypolish_options)
        step = Step(str(rule), polypolish, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_object(polypolish.informs, Path(output.INFORMS))

rule polishing_copy_fasta_polca:
    input:
        FASTA = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'polypolish' / 'polished.fasta'
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'polca' / 'input_genome.fasta'
    shell:
        """
        cp {input.FASTA} {output.FASTA}
        """

rule polishing_samtools_index_polca:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA_REF = rules.polishing_copy_fasta_polca.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'polca' / 'input_genome.fasta.fai'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'polca'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex(camel)
        samtools.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule), samtools, camel, Path(str(params.running_dir)))
        step.run_step()

rule polishing_polca:
    """
    Then polishing with Polca.
    """
    input:
        FQ_dict = Path(config['working_dir']) / 'fq_dict.io',
        FASTA = rules.polishing_copy_fasta_polca.output.FASTA,
        INDEX = rules.polishing_samtools_index_polca.output.INDEX_GENOME_PREFIX
    output:
        # FASTA = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'polca' / 'input_genome.fasta.PolcaCorrected.fa',
        FASTA = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'polca' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'polca' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'polca',
        polca_options = config.get('polishing', {}).get('polca', {})
    threads: 8
    run:
        from camel.app.tools.polca.polca import Polca
        polca = Polca(camel)
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict),'illumina')
        polca.add_input_files({'FASTQ_PE': fq_in.pe, 'FASTA': [ToolIOFile(Path(input.FASTA))]})
        polca.update_parameters(**params.polca_options, threads=threads)
        step = Step(str(rule), polca, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(polca, output)


rule hybrid_assembly_filter_contig_length:
    """
    Filters out the small contigs.
    """
    input:
        FASTA = rules.polishing_polca.output.FASTA
    output:
        FASTA = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'output' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / short_read_polishing.OUTPUT_ASSEMBLY_FILTERING_INFORMS
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'filtering',
        min_contig_length = config['assembly'].get('min_contig_length', 0) if 'assembly' in config else 0
    run:
        from camel.app.tools.seqtk.seqtkseq import SeqtkSeq
        seqtk = SeqtkSeq(camel)
        SnakemakeUtils.add_pickle_inputs(seqtk, input)
        step = Step(str(rule), seqtk, camel, Path(str(params.running_dir)))
        seqtk.update_parameters(output_filename='assembly_filtered.fasta', min_length=params.min_contig_length)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqtk, output)

rule hybrid_assembly_quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = rules.hybrid_assembly_filter_contig_length.output.FASTA
    output:
        TSV = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'quast' / 'tsv.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'quast'
    run:
        from camel.app.tools.quast.quast import Quast
        quast = Quast(camel)
        SnakemakeUtils.add_pickle_inputs(quast, input)
        step = Step(str(rule), quast, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast, output)

rule hybrid_assembly_quast_extract_informs:
    """
    Extracts the information from the QUAST output file.
    """
    input:
        TSV = rules.hybrid_assembly_quast.output.TSV
    output:
        INFORMS = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'quast' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'quast'
    run:
        from camel.app.tools.quast.quastinformextractor import QuastInformExtractor
        quast_inform_extractor = QuastInformExtractor(camel)
        SnakemakeUtils.add_pickle_inputs(quast_inform_extractor, input)
        step = Step(str(rule), quast_inform_extractor, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast_inform_extractor, output)

rule hybrid_assembly_report:
    """
    Creates the HTML report for the assembly.
    """
    input:
        FASTA_Raw = rules.polishing_polca.output.FASTA,
        FASTA_Contig = rules.hybrid_assembly_filter_contig_length.output.FASTA,
        INFORMS_spades = rules.polishing_polca.output.INFORMS,
        INFORMS_quast = rules.hybrid_assembly_quast_extract_informs.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / short_read_polishing.OUTPUT_ASSEMBLY_REPORT
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'report',
        sample_name = config['sample_name']
    run:
        from camel.app.tools.pipelines.assembly.htmlreporterassembly import HtmlReporterAssembly
        from camel.app.io.tooliovalue import ToolIOValue
        reporter = HtmlReporterAssembly(camel)
        reporter.add_input_files({'SAMPLE_NAME': [ToolIOValue(params.sample_name)],
                                  'ASSEMBLER': [ToolIOValue('Hybrid assembly')]})
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule hybrid_assembly_dump_summary_info:
    """
    Dumps the summary information from the assembly pipeline.
    """
    input:
        INFORMS_quast = rules.hybrid_assembly_quast_extract_informs.output.INFORMS
    output:
        Path(config['working_dir']) / short_read_polishing.OUTPUT_ASSEMBLY_SUMMARY
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'summary'
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

rule hybrid_assembly_bt2_index:
    """
    Creates a bowtie2 index for the assembly.
    """
    input:
        FASTA_REF = rules.hybrid_assembly_filter_contig_length.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'bowtie2' / 'genome_prefix.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'bowtie2'
    run:
        from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index
        bowtie2_index = Bowtie2Index(camel)
        step = Step(str(rule), bowtie2_index, camel, Path(str(params.running_dir)))
        SnakemakeUtils.add_pickle_inputs(bowtie2_index, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bowtie2_index, output)

rule hybrid_assembly_bt2_map:
    """
    Maps the reads against the assembled contigs.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
        INDEX_GENOME_PREFIX = rules.hybrid_assembly_bt2_index.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'bowtie2' / 'sam.io',
        INFORMS = Path(config['working_dir']) / short_read_polishing.OUTPUT_ASSEMBLY_MAPPING_INFORMS
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'bowtie2',
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        bowtie2_map = Bowtie2Map(camel)
        step = Step(str(rule), bowtie2_map, camel, Path(str(params.running_dir)))
        bowtie2_map.add_input_files(SnakePipelineUtils.extracts_fq_input(
            Path(input.IO), key_se='FASTQ_SE', drop_empty=True, read_type=params.read_type))
        SnakemakeUtils.add_pickle_input(bowtie2_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step.run_step()
        bowtie2_map.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_tool_outputs(bowtie2_map, output)

rule hybrid_assembly_bt2_sam_to_bam:
    """
    Converts the SAM file generated by bowtie2 to BAM format.
    """
    input:
        SAM = rules.hybrid_assembly_bt2_map.output.SAM
    output:
        BAM = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'bowtie2' / 'bam.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'bowtie2'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        samtools_view = SamtoolsView(camel)
        step = Step(str(rule), samtools_view, camel, Path(str(params.running_dir)))
        SnakemakeUtils.add_pickle_inputs(samtools_view, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_view, output)

rule hybrid_assembly_bt2_sort_bam:
    """
    Sorts the BAM alignment.
    """
    input:
        BAM = rules.hybrid_assembly_bt2_sam_to_bam.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'bowtie2' / 'bam_sorted.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'bowtie2'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort(camel)
        step = Step(str(rule), samtools_sort, camel, Path(str(params.running_dir)))
        SnakemakeUtils.add_pickle_inputs(samtools_sort, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_sort, output)

rule hybrid_assembly_samtools_depth:
    """
    Runs samtools depth on the BAM file of the reads mapped to the assembly.
    """
    input:
        BAM = rules.hybrid_assembly_bt2_sam_to_bam.output.BAM
    output:
        TSV = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'samtools_depth' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / short_read_polishing.OUTPUT_ASSEMBLY_DEPTH_INFORMS
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'samtools_depth'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth(camel)
        step = Step(str(rule), samtools_depth, camel, Path(str(params.running_dir)))
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        step.run_step()
        samtools_depth.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)

rule hybrid_assembly_nanopore_map_reads:
    """
    Maps the nanopore reads to the hybrid assembly.
    """
    input:
        FQ = Path(config['working_dir']) / 'fq_dict.io',
        FASTA = rules.hybrid_assembly_filter_contig_length.output.FASTA
    output:
        SAM = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'minimap2' / 'sam.io',
        INFORMS = Path(config['working_dir']) / short_read_polishing.OUTPUT_ASSEMBLY_NANOPORE_MAPPING_INFORMS
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'minimap2'
    run:
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        minimap2 = Minimap2Mapping(Camel.get_instance())
        SnakemakeUtils.add_pickle_input(minimap2, 'FASTA', Path(input.FASTA))
        minimap2.add_input_files(SnakePipelineUtils.extracts_fq_input(
            Path(input.FQ), key_se='FASTQ', read_type='SE'))
        step = Step(str(rule), minimap2, Camel.get_instance(), Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(minimap2, output)

rule hybrid_assembly_nanopore_sam_to_bam:
    """
    Converts the SAM file generated by minimap2 to BAM format.
    """
    input:
        SAM = rules.hybrid_assembly_nanopore_map_reads.output.SAM
    output:
        BAM = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'minimap2' / 'bam.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'minimap2'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        samtools_view = SamtoolsView(camel)
        step = Step(str(rule), samtools_view, camel, Path(str(params.running_dir)))
        SnakemakeUtils.add_pickle_inputs(samtools_view, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_view, output)

rule hybrid_assembly_nanopore_sort_bam:
    """
    Sorts the BAM alignment.
    """
    input:
        BAM = rules.hybrid_assembly_nanopore_sam_to_bam.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'minimap2' / 'bam_sorted.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'minimap2'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort(camel)
        step = Step(str(rule), samtools_sort, camel, Path(str(params.running_dir)))
        SnakemakeUtils.add_pickle_inputs(samtools_sort, input)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_sort, output)

rule hybrid_assembly_nanopore_samtools_flagstat:
    """
    Runs samtools flagstat to determine the mapping rate.
    """
    input:
        BAM = rules.hybrid_assembly_nanopore_sort_bam.output.BAM
    output:
        INFORMS = Path(config['working_dir']) / short_read_polishing.OUTPUT_ASSEMBLY_NANOPORE_MAPPING_RATE_INFORMS
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'samtools_flagstat'
    run:
        from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
        flagstat = SamtoolsFlagstat(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(flagstat, input)
        flagstat.run(Path(str(params.running_dir)).absolute())
        # Calculate mapping rate
        flagstat.informs['mapping_perc'] = 100 * flagstat.informs['mapped'][0] / flagstat.informs['total'][0]
        SnakemakeUtils.dump_tool_outputs(flagstat, output)

rule hybrid_assembly_nanopore_samtools_depth:
    """
    Runs samtools depth on the BAM file of the reads mapped to the assembly.
    """
    input:
        BAM = rules.hybrid_assembly_nanopore_sort_bam.output.BAM
    output:
        TSV = Path(config['working_dir']) / 'polishing' / '{assembly_type}' / 'samtools_depth_nanopore' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / short_read_polishing.OUTPUT_ASSEMBLY_NANOPORE_DEPTH_INFORMS
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'polishing' / wildcards.assembly_type / 'samtools_depth'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth(camel)
        step = Step(str(rule), samtools_depth, camel, Path(str(params.running_dir)))
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        step.run_step()
        samtools_depth.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)
