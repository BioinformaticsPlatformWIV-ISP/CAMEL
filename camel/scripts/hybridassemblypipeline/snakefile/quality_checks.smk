from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

camel = Camel.get_instance()

dictionary_references = {}

rule quast_final:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta'
    output:
        TSV = Path(config['working_dir']) / 'qc' / 'quast' / 'report.tsv'
    params:
        running_dir = Path(config['working_dir']) / 'qc' / 'quast'
    run:
        from camel.app.tools.quast.quast import Quast
        quast = Quast(camel)
        quast.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA))]})
        step = Step(str(rule), quast, camel, params.running_dir)
        step.run_step()

rule samtools_index_short_qc:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA_REF = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta'
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'polishing' / 'polypolish' /  'polished.fasta.fai'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polypolish'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex(camel)
        samtools.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule),samtools,camel,params.running_dir,config)
        step.run_step()

rule bwa_index_qc:
    """
    Creates a bwa index for the assembly.
    """
    input:
        FASTA_REF = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta'
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'polishing' / 'polypolish' /  'genome_prefix.io'
    params:
        running_dir = Path(config['working_dir']) / 'polishing' / 'polypolish'
    run:
        from camel.app.tools.bwa.bwaindex import BWAIndex
        bwa = BWAIndex(camel)
        bwa.add_input_files({'FASTA_REF': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule), bwa, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bwa, output)

rule read_mapping_qc:
    """
    Maps the reads against the assembly.
    """
    input:
        FQ_1P=Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_1P.fastq.gz",
        FQ_2P=Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_2P.fastq.gz",
        INDEX_GENOME_PREFIX_BWA = rules.bwa_index_qc.output.INDEX_GENOME_PREFIX,
        INDEX_GENOME_PREFIX_SAMTOOLS = rules.samtools_index_short_qc.output.INDEX_GENOME_PREFIX,
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta',
        INDEX_GENOME_PREFIX = rules.bwa_index_qc.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'qc' / 'read_mapping' / 'bwa_readmap.sam'
    params:
        running_dir = Path(config['working_dir']) / 'qc' / 'read_mapping'
    threads: 4
    run:
        from camel.app.tools.bwa.bwamap import BWAMap
        bwa_map = BWAMap(camel)
        bwa_map.add_input_files({'FASTQ_PE': [ToolIOFile(Path(input.FQ_1P)), ToolIOFile(Path(input.FQ_2P))]})
        SnakemakeUtils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step = Step(str(rule), bwa_map, camel, params.running_dir)
        step.run_step()

rule read_mapping_qc_longreads:
    """
    Maps the reads against the assembly.
    """
    input:
        FQ=Path(config['working_dir']) / 'trimming' / 'ont' / '{}_SE.fastq.gz'.format(config['name']),
        INDEX_GENOME_PREFIX_BWA = rules.bwa_index_qc.output.INDEX_GENOME_PREFIX,
        INDEX_GENOME_PREFIX_SAMTOOLS = rules.samtools_index_short_qc.output.INDEX_GENOME_PREFIX,
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta',
        INDEX_GENOME_PREFIX = rules.bwa_index_qc.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'qc' / 'read_mapping' / 'minimap2_readmap.sam'
    params:
        running_dir = Path(config['working_dir']) / 'qc' / 'read_mapping'
    threads: 4
    run:
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        minimap2 = Minimap2Mapping(camel)
        minimap2.add_input_files({'FASTQ': [ToolIOFile(Path(input.FQ))], 'FASTA':[ToolIOFile(Path(input.FASTA))]})
        minimap2.update_parameters(output_filename='minimap2_readmap.sam')
        step = Step(str(rule), minimap2, camel, params.running_dir, config)
        step.run_step()

rule sam_to_bam_qc:
    """
    Converts SAM to BAM.
    """
    input:
        SAM = rules.read_mapping_qc.output.SAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / 'read_mapping' / 'bwa_readmap.bam'
    params:
        running_dir = Path(config['working_dir']) / 'qc' / 'read_mapping'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        samtools_view = SamtoolsView(camel)
        samtools_view.add_input_files({'SAM': [ToolIOFile(Path(input.SAM))]})
        samtools_view.update_parameters(output_filename='bwa_readmap.bam')
        step = Step(rule, samtools_view, camel, params.running_dir, config)
        step.run_step()

rule bam_sorting_qc:
    """
    Sorts the BAM alignment.
    """
    input:
        BAM = rules.sam_to_bam_qc.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / 'read_mapping' / 'bwa_readmap.sorted.bam'
    params:
        running_dir = Path(config['working_dir']) / 'qc' / 'read_mapping'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort(camel)
        samtools_sort.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        samtools_sort.update_parameters(output_filename='bwa_readmap.sorted.bam')
        step = Step(rule, samtools_sort, camel, params.running_dir, config)
        step.run_step()

rule bam_indexing_qc:
    """
    Index the bam file.
    """
    input:
        BAM = rules.bam_sorting_qc.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / 'read_mapping' / 'bwa-index.io'
    params:
        running_dir = Path(config['working_dir']) /  'qc' / 'read_mapping'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        samtools_index = SamtoolsIndex(camel)
        samtools_index.add_input_files({'BAM':[ToolIOFile(Path(input.BAM))]})
        step = Step(str(rule), samtools_index, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_index, output)

rule sam_to_bam_qc_longreads:
    """
    Converts SAM to BAM.
    """
    input:
        SAM = rules.read_mapping_qc_longreads.output.SAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / 'read_mapping' / 'minimap2_readmap.bam'
    params:
        running_dir = Path(config['working_dir']) / 'qc' / 'read_mapping'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        samtools_view = SamtoolsView(camel)
        samtools_view.add_input_files({'SAM': [ToolIOFile(Path(input.SAM))]})
        samtools_view.update_parameters(output_filename='minimap2_readmap.bam')
        step = Step(rule, samtools_view, camel, params.running_dir, config)
        step.run_step()

rule bam_sorting_qc_longreads:
    """
    Sorts the BAM alignment.
    """
    input:
        BAM = rules.sam_to_bam_qc.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / 'read_mapping' / 'minimap2_readmap.sorted.bam'
    params:
        running_dir = Path(config['working_dir']) / 'qc' / 'read_mapping'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort(camel)
        samtools_sort.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        samtools_sort.update_parameters(output_filename='minimap2_readmap.sorted.bam')
        step = Step(rule, samtools_sort, camel, params.running_dir, config)
        step.run_step()

rule bam_indexing_qc_longreads:
    """
    Index the bam file.
    """
    input:
        BAM = rules.bam_sorting_qc_longreads.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / 'read_mapping' / 'bwa-index-longreads.io'
    params:
        running_dir = Path(config['working_dir']) /  'qc' / 'read_mapping'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        samtools_index = SamtoolsIndex(camel)
        samtools_index.add_input_files({'BAM':[ToolIOFile(Path(input.BAM))]})
        step = Step(str(rule), samtools_index, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_index, output)

rule freebayes_qc:
    """
    Checks for small variants in the final polished assembly.
    """
    input:
        BAM = rules.bam_sorting_qc.output.BAM,
        BAM_INDEX = rules.bam_indexing_qc.output.BAM,
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta',
        FASTA_INDEX = rules.samtools_index_short_qc.output.INDEX_GENOME_PREFIX
    output:
        VCF =  Path(config['working_dir']) / 'qc' / 'freebayes' / 'variants.vcf'
    params:
        running_dir = Path(config['working_dir']) /  'qc' / 'freebayes'
    run:
        from camel.app.tools.freebayes.freebayes import Freebayes
        freebayes = Freebayes(camel)
        freebayes.add_input_files({'BAM':[ToolIOFile(Path(input.BAM))], 'FASTA':[ToolIOFile(Path(input.FASTA))]})
        step = Step(rule, freebayes, camel, params.running_dir, config)
        step.run_step()

rule sniffles_qc:
    """
    Checks for structural variants in the final polished assembly using long reads.
    """
    input:
        BAM = rules.bam_sorting_qc_longreads.output.BAM,
        BAM_INDEX = rules.bam_indexing_qc_longreads.output.BAM,
        FASTA= Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta',
        FASTA_INDEX= rules.samtools_index_short_qc.output.INDEX_GENOME_PREFIX
    output:
        VCF = Path(config['working_dir']) / 'qc' / 'sniffles' / 'variants.vcf'
    params:
        running_dir = Path(config['working_dir']) /  'qc' / 'sniffles'
    run:
        from camel.app.tools.sniffles.sniffles import Sniffles
        sniffles = Sniffles(camel)
        sniffles.add_input_files({'BAM':[ToolIOFile(Path(input.BAM))], 'FASTA':[ToolIOFile(Path(input.FASTA))]})
        step = Step(rule, sniffles, camel, params.running_dir, config)
        step.run_step()

rule clair3_qc:
    """
    Checks for small variants in the final polished assembly using long reads.
    """
    input:
        BAM = rules.bam_sorting_qc_longreads.output.BAM,
        BAM_INDEX = rules.bam_indexing_qc_longreads.output.BAM,
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta',
        FASTA_INDEX = rules.samtools_index_short_qc.output.INDEX_GENOME_PREFIX
    output:
        VCF = Path(config['working_dir']) / 'qc' / 'clair3_output' / 'merge_output.vcf.gz'
    params:
        running_dir = Path(config['working_dir']) /  'qc'
    run:
        from camel.app.tools.clair3.clair3 import Clair3
        clair3 = Clair3(camel)
        clair3.add_input_files({'BAM':[ToolIOFile(Path(input.BAM))], 'FASTA':[ToolIOFile(Path(input.FASTA))]})
        step = Step(rule, clair3, camel, params.running_dir, config)
        step.run_step()

rule ale:
    """
    Generates ALE QC report for the final assembly
    """
    input:
        SAM = rules.read_mapping_qc.output.SAM,
        FASTA = Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta',
        FASTA_INDEX = rules.samtools_index_short_qc.output.INDEX_GENOME_PREFIX
    output:
        ALE = Path(config['working_dir']) / 'qc' / 'ale_illumina' / 'ALE.ale'
    params:
        running_dir = Path(config['working_dir']) / 'qc' / 'ale_illumina'
    run:
        from camel.app.tools.ale.ale import ALE
        ale = ALE(camel)
        ale.add_input_files({'SAM':[ToolIOFile(Path(input.SAM))], 'FASTA':[ToolIOFile(Path(input.FASTA))]})
        step = Step(rule, ale, camel, params.running_dir, config)
        step.run_step()

rule ale2wiggle:
    """
    Generates wiggle files from ALE for IGV visualization.
    """
    input:
        ALE = rules.ale.output.ALE
    output:
        TSV_1 = Path(config['working_dir']) / 'qc' / 'ale_illumina' / 'ALE.ale-depth.wig',
        TSV_2 = Path(config['working_dir']) / 'qc' / 'ale_illumina' / 'ALE.ale-kmer.wig',
        TSV_3 = Path(config['working_dir']) / 'qc' / 'ale_illumina' / 'ALE.ale-insert.wig',
        TSV_4 = Path(config['working_dir']) / 'qc' / 'ale_illumina' / 'ALE.ale-place.wig'
    params:
        running_dir = Path(config['working_dir']) / 'qc' / 'ale_illumina'
    run:
        from camel.app.tools.ale.ale2wiggle import ALE2Wiggle
        ale = ALE2Wiggle(camel)
        ale.add_input_files({'ALE':[ToolIOFile(Path(input.ALE))]})
        step = Step(rule, ale, camel, params.running_dir, config)
        step.run_step()
