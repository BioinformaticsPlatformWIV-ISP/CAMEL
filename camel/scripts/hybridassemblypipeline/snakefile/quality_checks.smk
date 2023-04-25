import pickle
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

camel = Camel.get_instance()

dictionary_references = {'medaka':[Path(config['working_dir']) / 'medaka' / 'consensus.fasta'],
                         'polca':[Path(config['working_dir']) / 'polishing' / 'polca' / 'consensus.fasta.PolcaCorrected.fa'],
                         'polypolish':[Path(config['working_dir']) / 'polishing' / 'polypolish' / 'polished.fasta'],
                         'unicycler':[Path(config['working_dir']) / 'unicycler' / 'assembly.fasta']}

rule check_qc:
    """
    Checks that the qc files are generated for each step of the assembly.
    """
    input:
        FASTA = [Path(config['working_dir']) / 'qc' / f'{name}' / 'ale_illumina' / 'ALE.ale-depth.wig' for name in dictionary_references.keys()]
    output:
        Path(config['working_dir'] / config['output'])
    shell:
        """
        touch {output}
        """

rule copy_fasta_file:
    """
    Moves the fasta file to qc location.
    """
    input:
        FASTA = lambda wildcards: dictionary_references[wildcards.name]
    output:
        FASTA = Path(config['working_dir']) / 'qc' / '{name}' / 'consensus.fasta'
    shell:
        """
        cp {input.FASTA} {output.FASTA}
        """

rule quast_final:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = rules.copy_fasta_file.output.FASTA
    output:
        TSV = Path(config['working_dir']) / 'qc' / '{name}' / 'quast' / 'report.tsv'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'quast'
    run:
        from camel.app.tools.quast.quast import Quast
        dir_working = Path(str(params.running_dir)).absolute()
        quast = Quast(camel)
        quast.add_input_files({'FASTA': [ToolIOFile(Path(str(input.FASTA)))]})
        step = Step(str(rule), quast, camel, dir_working)
        step.run_step()

rule parse_quast_output:
    """
    Parses the quast output into a json file.
    """
    input:
        TSV = rules.quast_final.output.TSV
    output:
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / 'quast' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'quast'
    run:
        from camel.app.tools.quast.quastinformextractor import QuastInformExtractor
        dir_working = Path(str(params.running_dir)).absolute()
        quast_inform_extractor = QuastInformExtractor(camel)
        quast_inform_extractor.add_input_files({'TSV': [ToolIOFile(Path(input.TSV))]})
        step = Step(str(rule), quast_inform_extractor, camel, dir_working, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast_inform_extractor, output)

rule samtools_index_short_qc:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA_REF = rules.copy_fasta_file.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'qc' / '{name}' / 'consensus.fasta.fai'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        dir_working = Path(str(params.running_dir)).absolute()
        samtools = SamtoolsFastaIndex(camel)
        samtools.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule),samtools,camel,dir_working,config)
        step.run_step()

rule bwa_index_qc:
    """
    Creates a bwa index for the assembly.
    """
    input:
        FASTA_REF = rules.copy_fasta_file.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'qc' / '{name}' /  'genome_prefix.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}'
    run:
        from camel.app.tools.bwa.bwaindex import BWAIndex
        dir_working = Path(str(params.running_dir)).absolute()
        bwa = BWAIndex(camel)
        bwa.add_input_files({'FASTA_REF': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule), bwa, camel, dir_working, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bwa, output)

rule read_mapping_qc:
    """
    Maps the reads against the assembly.
    """
    input:
        FQ_1P = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_1P.fastq.gz",
        FQ_2P = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_2P.fastq.gz",
        INDEX_GENOME_PREFIX_BWA = rules.bwa_index_qc.output.INDEX_GENOME_PREFIX,
        INDEX_GENOME_PREFIX_SAMTOOLS = rules.samtools_index_short_qc.output.INDEX_GENOME_PREFIX,
        FASTA = rules.copy_fasta_file.output.FASTA,
        INDEX_GENOME_PREFIX = rules.bwa_index_qc.output.INDEX_GENOME_PREFIX
    output:
        SAM =  Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'bwa_readmap.sam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping'
    threads: 4
    run:
        from camel.app.tools.bwa.bwamap import BWAMap
        dir_working = Path(str(params.running_dir)).absolute()
        bwa_map = BWAMap(camel)
        bwa_map.add_input_files({'FASTQ_PE': [ToolIOFile(Path(input.FQ_1P)), ToolIOFile(Path(input.FQ_2P))]})
        SnakemakeUtils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step = Step(str(rule), bwa_map, camel, dir_working)
        step.run_step()

rule read_mapping_qc_longreads:
    """
    Maps the reads against the assembly.
    """
    input:
        FQ = Path(config['working_dir']) / 'trimming' / 'ont' / '{}_SE.fastq.gz'.format(config['name']),
        INDEX_GENOME_PREFIX_BWA = rules.bwa_index_qc.output.INDEX_GENOME_PREFIX,
        INDEX_GENOME_PREFIX_SAMTOOLS = rules.samtools_index_short_qc.output.INDEX_GENOME_PREFIX,
        FASTA = rules.copy_fasta_file.output.FASTA,
        INDEX_GENOME_PREFIX = rules.bwa_index_qc.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'minimap2_readmap.sam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping'
    threads: 4
    run:
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        dir_working = Path(str(params.running_dir)).absolute()
        minimap2 = Minimap2Mapping(camel)
        minimap2.add_input_files({'FASTQ': [ToolIOFile(Path(input.FQ))], 'FASTA':[ToolIOFile(Path(input.FASTA))]})
        minimap2.update_parameters(output_filename='minimap2_readmap.sam')
        step = Step(str(rule), minimap2, camel, dir_working, config)
        step.run_step()

rule sam_to_bam_qc:
    """
    Converts SAM to BAM.
    """
    input:
        SAM = rules.read_mapping_qc.output.SAM
    output:
        BAM =  Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'bwa_readmap.bam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_view = SamtoolsView(camel)
        samtools_view.add_input_files({'SAM': [ToolIOFile(Path(input.SAM))]})
        samtools_view.update_parameters(output_filename='bwa_readmap.bam')
        step = Step(str(rule), samtools_view, camel, dir_working, config)
        step.run_step()

rule bam_sorting_qc:
    """
    Sorts the BAM alignment.
    """
    input:
        BAM = rules.sam_to_bam_qc.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'bwa_readmap.sorted.bam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_sort = SamtoolsSort(camel)
        samtools_sort.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        samtools_sort.update_parameters(output_filename='bwa_readmap.sorted.bam')
        step = Step(str(rule), samtools_sort, camel, dir_working, config)
        step.run_step()

rule bam_indexing_qc:
    """
    Index the bam file.
    """
    input:
        BAM = rules.bam_sorting_qc.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'bwa-index.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_index = SamtoolsIndex(camel)
        samtools_index.add_input_files({'BAM':[ToolIOFile(Path(input.BAM))]})
        step = Step(str(rule), samtools_index, camel, dir_working, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_index, output)

rule sam_to_bam_qc_longreads:
    """
    Converts SAM to BAM.
    """
    input:
        SAM = rules.read_mapping_qc_longreads.output.SAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'minimap2_readmap.bam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_view = SamtoolsView(camel)
        samtools_view.add_input_files({'SAM': [ToolIOFile(Path(input.SAM))]})
        samtools_view.update_parameters(output_filename='minimap2_readmap.bam')
        step = Step(str(rule), samtools_view, camel, dir_working, config)
        step.run_step()

rule bam_sorting_qc_longreads:
    """
    Sorts the BAM alignment.
    """
    input:
        BAM = rules.sam_to_bam_qc.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'minimap2_readmap.sorted.bam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_sort = SamtoolsSort(camel)
        samtools_sort.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        samtools_sort.update_parameters(output_filename='minimap2_readmap.sorted.bam')
        step = Step(str(rule), samtools_sort, camel, dir_working, config)
        step.run_step()

rule bam_indexing_qc_longreads:
    """
    Index the bam file.
    """
    input:
        BAM = rules.bam_sorting_qc_longreads.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'bwa-index-longreads.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_index = SamtoolsIndex(camel)
        samtools_index.add_input_files({'BAM':[ToolIOFile(Path(input.BAM))]})
        step = Step(str(rule), samtools_index, camel, dir_working, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_index, output)

rule freebayes_qc:
    """
    Checks for small variants in the final polished assembly.
    """
    input:
        BAM = rules.bam_sorting_qc.output.BAM,
        BAM_INDEX = rules.bam_indexing_qc.output.BAM,
        FASTA = rules.copy_fasta_file.output.FASTA,
        FASTA_INDEX = rules.samtools_index_short_qc.output.INDEX_GENOME_PREFIX
    output:
        VCF =  Path(config['working_dir']) / 'qc' / '{name}' / 'freebayes' / 'variants.vcf'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'freebayes',
        freebayes_options = config.get('freebayes',{})
    run:
        from camel.app.tools.freebayes.freebayes import Freebayes
        dir_working = Path(str(params.running_dir)).absolute()
        freebayes = Freebayes(camel)
        freebayes.add_input_files({'BAM':[ToolIOFile(Path(input.BAM))], 'FASTA':[ToolIOFile(Path(input.FASTA))]})
        freebayes.update_parameters(**params.freebayes_options)
        step = Step(str(rule), freebayes, camel, dir_working, config)
        step.run_step()

rule parse_freebayes_vcf:
    """
    Counts the number of variants in the freebayes vcf.
    """
    input:
        VCF = rules.freebayes_qc.output.VCF
    output:
        TSV = Path(config['working_dir']) / 'qc' / '{name}' / 'freebayes' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'freebayes'
    run:
        from camel.app.components.vcf.vcfutils import VCFUtils
        summary_json = {'number_of_variants': VCFUtils.count_variants(Path(input.VCF)),
                        'type_of_variants': VCFUtils.retrieve_variants(Path(input.VCF))}
        with open(output.TSV,'wb') as handle:
            pickle.dump(summary_json,handle)

rule sniffles_qc:
    """
    Checks for structural variants in the final polished assembly using long reads.
    """
    input:
        BAM = rules.bam_sorting_qc_longreads.output.BAM,
        BAM_INDEX = rules.bam_indexing_qc_longreads.output.BAM,
        FASTA= rules.copy_fasta_file.output.FASTA,
        FASTA_INDEX= rules.samtools_index_short_qc.output.INDEX_GENOME_PREFIX
    output:
        VCF = Path(config['working_dir']) / 'qc' / '{name}' / 'sniffles' / 'variants.vcf'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'sniffles'
    run:
        from camel.app.tools.sniffles.sniffles import Sniffles
        dir_working = Path(str(params.running_dir)).absolute()
        sniffles = Sniffles(camel)
        sniffles.add_input_files({'BAM':[ToolIOFile(Path(input.BAM))], 'FASTA':[ToolIOFile(Path(input.FASTA))]})
        step = Step(str(rule), sniffles, camel, dir_working, config)
        step.run_step()

rule clair3_qc:
    """
    Checks for small variants in the final polished assembly using long reads.
    """
    input:
        BAM = rules.bam_sorting_qc_longreads.output.BAM,
        BAM_INDEX = rules.bam_indexing_qc_longreads.output.BAM,
        FASTA = rules.copy_fasta_file.output.FASTA,
        FASTA_INDEX = rules.samtools_index_short_qc.output.INDEX_GENOME_PREFIX
    output:
        VCF = Path(config['working_dir']) / 'qc' / '{name}' / 'clair3_output' / 'merge_output.vcf.gz'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}',
        clair3_options= config.get('clair3',{})
    run:
        from camel.app.tools.clair3.clair3 import Clair3
        dir_working = Path(str(params.running_dir)).absolute()
        clair3 = Clair3(camel)
        clair3.add_input_files({'BAM':[ToolIOFile(Path(input.BAM))], 'FASTA':[ToolIOFile(Path(input.FASTA))]})
        clair3.update_parameters(**params.clair3_options)
        step = Step(str(rule), clair3, camel, dir_working, config)
        step.run_step()

rule parse_clair3_vcf:
    """
    Counts the number of variants in the clair3 vcf.
    """
    input:
        VCF = rules.clair3_qc.output.VCF
    output:
        TSV = Path(config['working_dir']) / 'qc' / '{name}' / 'clair3_output' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'clair3_output'
    run:
        from camel.app.components.vcf.vcfutils import VCFUtils
        summary_json = {'number_of_variants': VCFUtils.count_variants(Path(input.VCF)),
                        'type_of_variants': VCFUtils.retrieve_variants(Path(input.VCF))}
        with open(output.TSV,'wb') as handle:
            pickle.dump(summary_json,handle)

rule parse_sniffles_vcf:
    """
    Counts the number of variants in the sniffles vcf.
    """
    input:
        VCF = rules.sniffles_qc.output.VCF
    output:
        TSV = Path(config['working_dir']) / 'qc' / '{name}' / 'sniffles' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'sniffles'
    run:
        from camel.app.components.vcf.vcfutils import VCFUtils
        summary_json = {'number_of_variants': VCFUtils.count_variants(Path(input.VCF)),
                        'type_of_variants': VCFUtils.retrieve_variants(Path(input.VCF))}
        with open(output.TSV,'wb') as handle:
            pickle.dump(summary_json,handle)

rule ale:
    """
    Generates ALE QC report for the final assembly
    """
    input:
        SAM = rules.read_mapping_qc.output.SAM,
        FASTA = rules.copy_fasta_file.output.FASTA,
        FASTA_INDEX = rules.samtools_index_short_qc.output.INDEX_GENOME_PREFIX
    output:
        ALE = Path(config['working_dir']) / 'qc' / '{name}' / 'ale_illumina' / 'ALE.ale'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'ale_illumina'
    run:
        from camel.app.tools.ale.ale import ALE
        dir_working = Path(str(params.running_dir)).absolute()
        ale_report = ALE(camel)
        ale_report.add_input_files({'SAM':[ToolIOFile(Path(input.SAM))], 'FASTA':[ToolIOFile(Path(input.FASTA))]})
        step = Step(str(rule), ale_report, camel, dir_working, config)
        step.run_step()

rule ale2wiggle:
    """
    Generates wiggle files from ALE for IGV visualization.
    """
    input:
        ALE = rules.ale.output.ALE
    output:
        TSV_1 = Path(config['working_dir']) / 'qc' / '{name}' / 'ale_illumina' / 'ALE.ale-depth.wig',
        TSV_2 = Path(config['working_dir']) / 'qc' / '{name}' / 'ale_illumina' / 'ALE.ale-kmer.wig',
        TSV_3 = Path(config['working_dir']) / 'qc' / '{name}' / 'ale_illumina' / 'ALE.ale-insert.wig',
        TSV_4 = Path(config['working_dir']) / 'qc' / '{name}' / 'ale_illumina' / 'ALE.ale-place.wig'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'ale_illumina'
    run:
        from camel.app.tools.ale.ale2wiggle import ALE2Wiggle
        dir_working = Path(str(params.running_dir)).absolute()
        ale2wiggle_report = ALE2Wiggle(camel)
        ale2wiggle_report.add_input_files({'ALE':[ToolIOFile(Path(input.ALE))]})
        step = Step(str(rule), ale2wiggle_report, camel, dir_working, config)
        step.run_step()
