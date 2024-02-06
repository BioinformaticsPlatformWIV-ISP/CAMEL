from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.hybridassemblypipeline.snakefile import qc_hybrid

camel = Camel.get_instance()

rule qc_hybrid_samtools_index:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA = lambda wildcards: str(Path(config['working_dir']) / qc_hybrid.consensus_by_tool[wildcards.name])
    output:
        FASTA = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'fasta-index.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        dir_working = Path(str(params.running_dir)).absolute()
        samtools = SamtoolsFastaIndex(camel)
        SnakemakeUtils.add_pickle_inputs(samtools, input)
        step = Step(str(rule), samtools, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools, output)

rule qc_hybrid_quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        TSV = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'quast' / 'report.tsv',
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'quast' / 'commands.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'quast'
    run:
        from camel.app.tools.quast.quast import Quast
        dir_working = Path(str(params.running_dir)).absolute()
        quast = Quast(camel)
        SnakemakeUtils.add_pickle_inputs(quast, input)
        step = Step(str(rule), quast, camel, dir_working)
        step.run_step()
        quast.informs['_tag'] = wildcards.name
        SnakemakeUtils.dump_object(quast.informs, Path(output.INFORMS))

rule qc_hybrid_parse_quast_output:
    """
    Parses the quast output into a pickle.
    """
    input:
        TSV = rules.qc_hybrid_quast.output.TSV
    output:
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'quast' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'quast'
    run:
        from camel.app.tools.quast.quastinformextractor import QuastInformExtractor
        dir_working = Path(str(params.running_dir)).absolute()
        quast_inform_extractor = QuastInformExtractor(camel)
        quast_inform_extractor.add_input_files({'TSV': [ToolIOFile(Path(input.TSV))]})
        step = Step(str(rule), quast_inform_extractor, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast_inform_extractor, output)

rule qc_hybrid_quast_all_assemblies:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = [Path(config['working_dir']) / 'qc_hybrid' / f'{name}' / 'fasta-index.io' for name in config['assembly_steps']]
    output:
        TSV = Path(config['working_dir']) / 'qc_hybrid' / 'quast_combined' / 'report.tsv',
        HTML = Path(config['working_dir']) / 'qc_hybrid' / 'quast_combined' / 'report.html',
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / 'quast_combined' / 'commands.io'
    params:
        running_dir = Path(config['working_dir']) / 'qc_hybrid'  / 'quast_combined'
    run:
        from camel.app.tools.quast.quast import Quast
        dir_working = Path(str(params.running_dir)).absolute()
        quast = Quast(camel)
        quast.add_input_files({'FASTA': [SnakemakeUtils.load_object(Path(str(i)))[0] for i in input.FASTA]})
        step = Step(str(rule), quast, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_object(quast.informs, Path(output.INFORMS))

rule qc_hybrid_bwa_index:
    """
    Creates a bwa index for the assembly.
    """
    input:
        FASTA_REF = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'qc_hybrid' / '{name}' /  'genome_prefix.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}'
    run:
        from camel.app.tools.bwa.bwaindex import BWAIndex
        dir_working = Path(str(params.running_dir)).absolute()
        bwa = BWAIndex(camel)
        SnakemakeUtils.add_pickle_inputs(bwa, input)
        step = Step(str(rule), bwa, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bwa, output)

rule qc_hybrid_read_mapping_illumina:
    """
    Maps the short reads against the assembly.
    """
    input:
        FQ_dict = Path(config['working_dir']) / 'fq_dict.io',
        INDEX_GENOME_PREFIX_BWA = rules.qc_hybrid_bwa_index.output.INDEX_GENOME_PREFIX,
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        SAM = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'illumina' / 'bwa_readmap_sam.io',
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'illumina' / 'commands.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'read_mapping' / 'illumina'
    threads: 8
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.bwa.bwamap import BWAMap
        dir_working = Path(str(params.running_dir)).absolute()
        bwa_map = BWAMap(camel)
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict), 'illumina')
        bwa_map.add_input_files({'FASTQ_PE': fq_in.pe})
        bwa_map.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX_BWA))
        step = Step(str(rule), bwa_map, camel, dir_working)
        step.run_step()
        bwa_map.informs['_tag'] = wildcards.name
        SnakemakeUtils.dump_tool_outputs(bwa_map, output)

rule qc_hybrid_read_mapping_ont:
    """
    Maps the long reads against the assembly.
    """
    input:
        FQ = Path(config['working_dir']) / 'fq_dict.io',
        INDEX_GENOME_PREFIX_BWA = rules.qc_hybrid_bwa_index.output.INDEX_GENOME_PREFIX,
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        SAM = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'ont' / 'ont_mapping.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'read_mapping' / 'ont'
    threads: 8
    run:
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        dir_working = Path(str(params.running_dir)).absolute()
        minimap2 = Minimap2Mapping(camel)
        SnakemakeUtils.add_pickle_input(minimap2, 'FASTA', Path(input.FASTA))
        minimap2.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.FQ), key_se='FASTQ', read_type='SE'))
        step = Step(str(rule), minimap2, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(minimap2, output)

rule qc_hybrid_sam_to_bam_illumina:
    """
    Converts SAM to BAM from the short reads mapping.
    """
    input:
        SAM = rules.qc_hybrid_read_mapping_illumina.output.SAM
    output:
        BAM =  Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'illumina' / 'bwa_readmap_bam.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'read_mapping' / 'illumina'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_view = SamtoolsView(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_view, input)
        step = Step(str(rule), samtools_view, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_view, output)

rule qc_hybrid_bam_sorting_illumina:
    """
    Sorts the BAM alignment from the short reads mapping.
    """
    input:
        BAM = rules.qc_hybrid_sam_to_bam_illumina.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'illumina' / 'bwa_readmap_sorted.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'read_mapping' / 'illumina'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_sort = SamtoolsSort(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_sort, input)
        step = Step(str(rule), samtools_sort, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_sort, output)

rule qc_hybrid_bam_indexing_illumina:
    """
    Indexes the BAM file from the short reads mapping.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_illumina.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'illumina' / 'samtools-index.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'read_mapping' / 'illumina'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_index = SamtoolsIndex(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_index, input)
        step = Step(str(rule), samtools_index, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_index, output)

rule qc_hybrid_mapping_stats_illumina:
    """
    Retrieves mapping stats from the short read mappings.
    """
    input:
        # BAM = rules.qc_hybrid_bam_sorting_illumina.output.BAM,
        BAM = rules.qc_hybrid_bam_indexing_illumina.output.BAM
    output:
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'illumina' / 'flagstat.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'read_mapping' / 'illumina'
    run:
        from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_flagstat = SamtoolsFlagstat(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_flagstat, input)
        step = Step(str(rule), samtools_flagstat, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_flagstat, output)

rule qc_hybrid_samtools_depth_illumina:
    """
    Runs samtools depth on the BAM file from the short reads mapping.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_illumina.output.BAM,
        INDEX = rules.qc_hybrid_bam_indexing_illumina.output.BAM
    output:
        TSV = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'illumina' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'illumina' / 'samtools-depth.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'read_mapping' / 'illumina'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_depth = SamtoolsDepth(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        step = Step(str(rule), samtools_depth, camel, dir_working)
        step.run_step()
        samtools_depth.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)

rule qc_hybrid_sam_to_bam_ont:
    """
    Converts SAM to BAM from the long reads mapping.
    """
    input:
        SAM = rules.qc_hybrid_read_mapping_ont.output.SAM
    output:
        BAM = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'ont' / 'minimap2_bam.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'read_mapping' / 'ont'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_view = SamtoolsView(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_view, input)
        step = Step(str(rule), samtools_view, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_view, output)

rule qc_hybrid_bam_sorting_ont:
    """
    Sorts the BAM alignment from the long reads mapping.
    """
    input:
        BAM = rules.qc_hybrid_sam_to_bam_ont.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'ont' / 'minimap2_sorted.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'read_mapping' / 'ont'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_sort = SamtoolsSort(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_sort, input)
        step = Step(str(rule), samtools_sort, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_sort, output)

rule qc_hybrid_bam_indexing_ont:
    """
    Indexes the BAM file from the long reads mapping.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_ont.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'ont' / 'bwa-index-longreads.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'read_mapping' / 'ont'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_index = SamtoolsIndex(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_index, input)
        step = Step(str(rule), samtools_index, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_index, output)

rule qc_hybrid_mapping_stats_ont:
    """
    Retrieves mapping stats from the long read mapping.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_ont.output.BAM,
        INDEX = rules.qc_hybrid_bam_indexing_ont.output.BAM
    output:
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'ont' / 'flagstat-longreads.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'read_mapping' / 'ont'
    run:
        from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_flagstat = SamtoolsFlagstat(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_flagstat, input)
        step = Step(str(rule), samtools_flagstat, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_flagstat, output)

rule qc_hybrid_samtools_depth_ont:
    """
    Runs samtools depth on the BAM file of the long reads mapping.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_ont.output.BAM,
        INDEX = rules.qc_hybrid_bam_indexing_ont.output.BAM
    output:
        TSV = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'ont' / 'tsv-long.io',
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'read_mapping' / 'ont' / 'samtools-depth-long.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'read_mapping' / 'ont'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_depth = SamtoolsDepth(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        step = Step(str(rule), samtools_depth, camel, dir_working)
        step.run_step()
        samtools_depth.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)

rule qc_hybrid_freebayes:
    """
    Checks for small variants in the assembly.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_illumina.output.BAM,
        BAM_INDEX = rules.qc_hybrid_bam_indexing_illumina.output.BAM,
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        VCF =  Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'freebayes' / 'vcf.io',
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'freebayes' / 'commands.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'freebayes',
        freebayes_options = config.get('freebayes',{})
    run:
        from camel.app.tools.freebayes.freebayes import Freebayes
        dir_working = Path(str(params.running_dir)).absolute()
        freebayes = Freebayes(camel)
        SnakemakeUtils.add_pickle_inputs(freebayes, input)
        freebayes.update_parameters(**params.freebayes_options)
        step = Step(str(rule), freebayes, camel, dir_working)
        step.run_step()
        freebayes.informs['_tag'] = wildcards.name
        SnakemakeUtils.dump_tool_outputs(freebayes, output)

rule qc_hybrid_sniffles:
    """
    Checks for structural variants in the assembly using long reads.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_ont.output.BAM,
        BAM_INDEX = rules.qc_hybrid_bam_indexing_ont.output.BAM,
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        VCF = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'sniffles' / 'vcf.io',
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'sniffles' / 'commands.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'sniffles'
    threads: 4
    run:
        from camel.app.tools.sniffles.sniffles import Sniffles
        dir_working = Path(str(params.running_dir)).absolute()
        sniffles = Sniffles(camel)
        SnakemakeUtils.add_pickle_inputs(sniffles, input)
        sniffles.update_parameters(threads=threads)
        step = Step(str(rule), sniffles, camel, dir_working)
        step.run_step()
        sniffles.informs['_tag'] = wildcards.name
        SnakemakeUtils.dump_tool_outputs(sniffles, output)

rule qc_hybrid_clair3:
    """
    Checks for small variants in the assembly using long reads.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_ont.output.BAM,
        BAM_INDEX = rules.qc_hybrid_bam_indexing_ont.output.BAM,
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        VCF = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'clair3_output' / 'gzipped_vcf.io',
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'clair3_output' / 'commands.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}',
        clair3_options = config.get('clair3', {})
    threads: 8
    run:
        from camel.app.tools.clair3.clair3 import Clair3
        dir_working = Path(str(params.running_dir)).absolute()
        clair3 = Clair3(camel)
        SnakemakeUtils.add_pickle_inputs(clair3, input)
        clair3.update_parameters(
            **params.clair3_options, chunk_size=100_000, platform='ont', no_phasing=True, include_ctgs=True,
            threads=threads)
        step = Step(str(rule), clair3, camel, dir_working)
        step.run_step()
        clair3.informs['_tag'] = wildcards.name
        SnakemakeUtils.dump_tool_outputs(clair3, output)

rule qc_hybrid_unzip_clair3_vcf:
    """
    Unzips the output from clair3 to be compatible with the qc_add_vcf_info_to_informs rule.
    """
    input:
        VCF_GZ = rules.qc_hybrid_clair3.output.VCF
    output:
        VCF = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'clair3_output' / 'vcf.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'clair3_output'
    run:
        from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
        dir_working = Path(str(params.running_dir)).absolute()
        bcftools_view = BcftoolsView(camel)
        SnakemakeUtils.add_pickle_inputs(bcftools_view, input)
        bcftools_view.update_parameters(output_format='VCF', compress_output=False, output_filename='variants.vcf')
        step = Step(str(rule), bcftools_view, camel, dir_working)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(bcftools_view, output)

rule qc_hybrid_add_vcf_info_to_informs:
    """
    Parses the VCF file and adds the statistics to the informs.
    """
    input:
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / '{method}' / 'commands.io',
        VCF = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / '{method}' / 'vcf.io'
    output:
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / '{method}' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / wildcards.name / wildcards.method
    run:
        from camel.app.components.vcf.vcfutils import VCFUtils
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        if wildcards.method == 'sniffles':
            informs['nb_of_variants'] = sum(informs['variants'].values())
            informs['nb_of_snps'] = 0
            informs['nb_of_indels'] = informs['variants']['DEL'] + informs['variants']['INS']
            informs['nb_of_svs'] = informs['variants']['BND'] + informs['variants']['DUP'] + informs['variants']['INV']
        else:
            variants = VCFUtils.retrieve_variants(SnakemakeUtils.load_object(Path(input.VCF))[0].path)
            informs['nb_of_variants'] = len(variants)
            informs['nb_of_snps'] = sum(v.is_snp for v in variants)
            informs['nb_of_indels'] = sum(v.is_indel for v in variants)
            informs['nb_of_svs'] = 0
        SnakemakeUtils.dump_object(informs, Path(output.INFORMS))

rule qc_hybrid_ale:
    """
    Generates ALE QC report for the final assembly.
    """
    input:
        SAM = rules.qc_hybrid_read_mapping_illumina.output.SAM,
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        ALE = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'ale_illumina' / 'ALE.io',
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'ale_illumina' / 'informs-report.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'ale_illumina'
    run:
        from camel.app.tools.ale.ale import ALE
        dir_working = Path(str(params.running_dir)).absolute()
        ale_report = ALE(camel)
        SnakemakeUtils.add_pickle_inputs(ale_report, input)
        step = Step(str(rule), ale_report, camel, dir_working)
        step.run_step()
        ale_report.informs['_tag'] = f'{wildcards.name}'
        SnakemakeUtils.dump_tool_outputs(ale_report, output)

rule qc_hybrid_ale2wiggle:
    """
    Generates wiggle files from ALE for IGV visualization.
    """
    input:
        ALE = rules.qc_hybrid_ale.output.ALE
    output:
        TSV = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'ale_illumina' / 'wiggle.io',
        INFORMS = Path(config['working_dir']) / 'qc_hybrid' / '{name}' / 'ale_illumina' / 'commands-wiggle.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc_hybrid' / f'{wildcards.name}' / 'ale_illumina'
    run:
        from camel.app.tools.ale.ale2wiggle import ALE2Wiggle
        dir_working = Path(str(params.running_dir)).absolute()
        ale2wiggle_report = ALE2Wiggle(camel)
        SnakemakeUtils.add_pickle_inputs(ale2wiggle_report, input)
        step = Step(str(rule), ale2wiggle_report, camel, dir_working)
        step.run_step()
        ale2wiggle_report.informs['_tag'] = f'{wildcards.name}'
        SnakemakeUtils.dump_tool_outputs(ale2wiggle_report, output)
