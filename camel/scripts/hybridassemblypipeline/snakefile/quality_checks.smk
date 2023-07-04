from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.hybridassemblypipeline.snakefile import quality_checks

camel = Camel.get_instance()

rule qc_copy_fasta_file:
    """
    Moves the fasta file to the qc location.
    """
    input:
        FASTA = lambda wildcards: str(Path(config['working_dir']) / quality_checks.consensus_by_tool[wildcards.name])
    output:
        FASTA = Path(config['working_dir']) / 'qc' / '{name}' / 'consensus.fasta'
    shell:
        """
        cp {input.FASTA} {output.FASTA}
        """

rule qc_quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = rules.qc_copy_fasta_file.output.FASTA
    output:
        TSV = Path(config['working_dir']) / 'qc' / '{name}' / 'quast' / 'report.tsv',
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / 'quast' / 'commands.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'quast'
    run:
        from camel.app.tools.quast.quast import Quast
        dir_working = Path(str(params.running_dir)).absolute()
        quast = Quast(camel)
        quast.add_input_files({'FASTA': [ToolIOFile(Path(str(input.FASTA)))]})
        step = Step(str(rule), quast, camel, dir_working, config)
        step.run_step()
        quast.informs['_tag'] = wildcards.name
        SnakemakeUtils.dump_object(quast.informs, Path(output.INFORMS))

rule qc_parse_quast_output:
    """
    Parses the quast output into a pickle.
    """
    input:
        TSV = rules.qc_quast.output.TSV
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

rule qc_quast_all_assemblies:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = [Path(config['working_dir']) / 'qc' / f'{name}' / 'consensus.fasta' for name in quality_checks.consensus_by_tool.keys()]
    output:
        TSV = Path(config['working_dir']) / 'qc' / 'quast_combined' / 'report.tsv',
        HTML = Path(config['working_dir']) / 'qc' / 'quast_combined' / 'report.html',
        INFORMS = Path(config['working_dir']) / 'qc' / 'quast_combined' / 'commands.io'
    params:
        running_dir = Path(config['working_dir']) / 'qc'  / 'quast_combined'
    run:
        from camel.app.tools.quast.quast import Quast
        dir_working = Path(str(params.running_dir)).absolute()
        quast = Quast(camel)
        quast.add_input_files({'FASTA': [ToolIOFile(Path(str(i))) for i in input.FASTA]})
        step = Step(str(rule), quast, camel, dir_working, config)
        step.run_step()
        SnakemakeUtils.dump_object(quast.informs, Path(output.INFORMS))

rule qc_samtools_index:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA_REF = rules.qc_copy_fasta_file.output.FASTA
    output:
        INDEX_GENOME_PREFIX = Path(config['working_dir']) / 'qc' / '{name}' / 'consensus.fasta.fai'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        dir_working = Path(str(params.running_dir)).absolute()
        samtools = SamtoolsFastaIndex(camel)
        samtools.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA_REF))]})
        step = Step(str(rule), samtools, camel, dir_working, config)
        step.run_step()

rule qc_bwa_index:
    """
    Creates a bwa index for the assembly.
    """
    input:
        FASTA_REF = rules.qc_copy_fasta_file.output.FASTA
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

rule qc_read_mapping_illumina:
    """
    Maps the short reads against the assembly.
    """
    input:
        FQ_dict = Path(config['working_dir']) / 'trimming' / 'illumina' / 'fq_dict.io',
        INDEX_GENOME_PREFIX_BWA = rules.qc_bwa_index.output.INDEX_GENOME_PREFIX,
        INDEX_GENOME_PREFIX_SAMTOOLS = rules.qc_samtools_index.output.INDEX_GENOME_PREFIX,
        FASTA = rules.qc_copy_fasta_file.output.FASTA,
        INDEX_GENOME_PREFIX = rules.qc_bwa_index.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'illumina' / 'bwa_readmap.sam',
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'illumina' / 'commands.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping' / 'illumina'
    threads: 8
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.bwa.bwamap import BWAMap
        dir_working = Path(str(params.running_dir)).absolute()
        bwa_map = BWAMap(camel)
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict), 'illumina')
        bwa_map.add_input_files({'FASTQ_PE': fq_in.pe})
        bwa_map.update_parameters(threads=threads)
        SnakemakeUtils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step = Step(str(rule), bwa_map, camel, dir_working, config)
        step.run_step()
        bwa_map.informs['_tag'] = wildcards.name
        SnakemakeUtils.dump_object(bwa_map.informs, Path(output.INFORMS))

rule qc_read_mapping_ont:
    """
    Maps the long reads against the assembly.
    """
    input:
        FQ = Path(config['working_dir']) / 'trimming' / 'ont' / 'trimmed.fastq.gz',
        INDEX_GENOME_PREFIX_BWA = rules.qc_bwa_index.output.INDEX_GENOME_PREFIX,
        INDEX_GENOME_PREFIX_SAMTOOLS = rules.qc_samtools_index.output.INDEX_GENOME_PREFIX,
        FASTA = rules.qc_copy_fasta_file.output.FASTA,
        INDEX_GENOME_PREFIX = rules.qc_bwa_index.output.INDEX_GENOME_PREFIX
    output:
        SAM = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'ont' / 'minimap2_readmap.sam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping' / 'ont'
    threads: 8
    run:
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        dir_working = Path(str(params.running_dir)).absolute()
        minimap2 = Minimap2Mapping(camel)
        minimap2.add_input_files({'FASTQ': [ToolIOFile(Path(input.FQ))], 'FASTA': [ToolIOFile(Path(input.FASTA))]})
        minimap2.update_parameters(output_filename='minimap2_readmap.sam', threads=threads)
        step = Step(str(rule), minimap2, camel, dir_working, config)
        step.run_step()

rule qc_sam_to_bam_illumina:
    """
    Converts SAM to BAM from the short reads mapping.
    """
    input:
        SAM = rules.qc_read_mapping_illumina.output.SAM
    output:
        BAM =  Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'illumina' / 'bwa_readmap.bam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping' / 'illumina'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_view = SamtoolsView(camel)
        samtools_view.add_input_files({'SAM': [ToolIOFile(Path(input.SAM))]})
        samtools_view.update_parameters(output_filename='bwa_readmap.bam')
        step = Step(str(rule), samtools_view, camel, dir_working, config)
        step.run_step()

rule qc_bam_sorting_illumina:
    """
    Sorts the BAM alignment from the short reads mapping.
    """
    input:
        BAM = rules.qc_sam_to_bam_illumina.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'illumina' / 'bwa_readmap.sorted.bam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping' / 'illumina'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_sort = SamtoolsSort(camel)
        samtools_sort.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        samtools_sort.update_parameters(output_filename='bwa_readmap.sorted.bam')
        step = Step(str(rule), samtools_sort, camel, dir_working, config)
        step.run_step()

rule qc_bam_indexing_illumina:
    """
    Indexes the BAM file from the short reads mapping.
    """
    input:
        BAM = rules.qc_bam_sorting_illumina.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'illumina' / 'bwa-index.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping' / 'illumina'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_index = SamtoolsIndex(camel)
        samtools_index.add_input_files({'BAM':[ToolIOFile(Path(input.BAM))]})
        step = Step(str(rule), samtools_index, camel, dir_working, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_index, output)

rule qc_mapping_stats_illumina:
    """
    Retrieves mapping stats from the short read mappings.
    """
    input:
        BAM = rules.qc_bam_sorting_illumina.output.BAM,
        INDEX = rules.qc_bam_indexing_illumina.output.BAM
    output:
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'illumina' / 'flagstat.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping' / 'illumina'
    run:
        from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_flagstat = SamtoolsFlagstat(camel)
        samtools_flagstat.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        step = Step(str(rule), samtools_flagstat, camel, dir_working, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_flagstat, output)

rule qc_samtools_depth_illumina:
    """
    Runs samtools depth on the BAM file from the short reads mapping.
    """
    input:
        BAM = rules.qc_bam_sorting_illumina.output.BAM,
        INDEX = rules.qc_bam_indexing_illumina.output.BAM
    output:
        TSV = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'illumina' / 'tsv.io',
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'illumina' / 'samtools-depth.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping' / 'illumina'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_depth = SamtoolsDepth(camel)
        samtools_depth.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        step = Step(str(rule), samtools_depth, camel, dir_working, config)
        step.run_step()
        samtools_depth.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)

rule qc_sam_to_bam_ont:
    """
    Converts SAM to BAM from the long reads mapping.
    """
    input:
        SAM = rules.qc_read_mapping_ont.output.SAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'ont' / 'minimap2_readmap.bam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping' / 'ont'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_view = SamtoolsView(camel)
        samtools_view.add_input_files({'SAM': [ToolIOFile(Path(input.SAM))]})
        samtools_view.update_parameters(output_filename='minimap2_readmap.bam')
        step = Step(str(rule), samtools_view, camel, dir_working, config)
        step.run_step()

rule qc_bam_sorting_ont:
    """
    Sorts the BAM alignment from the long reads mapping.
    """
    input:
        BAM = rules.qc_sam_to_bam_ont.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'ont' / 'minimap2_readmap.sorted.bam'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping' / 'ont'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_sort = SamtoolsSort(camel)
        samtools_sort.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        samtools_sort.update_parameters(output_filename='minimap2_readmap.sorted.bam')
        step = Step(str(rule), samtools_sort, camel, dir_working, config)
        step.run_step()

rule qc_bam_indexing_ont:
    """
    Indexes the BAM file from the long reads mapping.
    """
    input:
        BAM = rules.qc_bam_sorting_ont.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'ont' / 'bwa-index-longreads.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping' / 'ont'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_index = SamtoolsIndex(camel)
        samtools_index.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        step = Step(str(rule), samtools_index, camel, dir_working, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_index, output)

rule qc_mapping_stats_ont:
    """
    Retrieves mapping stats from the long read mapping.
    """
    input:
        BAM = rules.qc_bam_sorting_ont.output.BAM,
        INDEX = rules.qc_bam_indexing_ont.output.BAM
    output:
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'ont' / 'flagstat-longreads.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping' / 'ont'
    run:
        from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_flagstat = SamtoolsFlagstat(camel)
        samtools_flagstat.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        step = Step(str(rule), samtools_flagstat, camel, dir_working, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_flagstat, output)

rule qc_samtools_depth_ont:
    """
    Runs samtools depth on the BAM file of the long reads mapping.
    """
    input:
        BAM = rules.qc_bam_sorting_ont.output.BAM,
        INDEX = rules.qc_bam_indexing_ont.output.BAM
    output:
        TSV = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'ont' / 'tsv-long.io',
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / 'read_mapping' / 'ont' / 'samtools-depth-long.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'read_mapping' / 'ont'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        dir_working = Path(str(params.running_dir)).absolute()
        samtools_depth = SamtoolsDepth(camel)
        samtools_depth.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        step = Step(str(rule), samtools_depth, camel, dir_working, config)
        step.run_step()
        samtools_depth.informs['_tag'] = 'Coverage calculation'
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)

rule qc_freebayes:
    """
    Checks for small variants in the assembly.
    """
    input:
        BAM = rules.qc_bam_sorting_illumina.output.BAM,
        BAM_INDEX = rules.qc_bam_indexing_illumina.output.BAM,
        FASTA = rules.qc_copy_fasta_file.output.FASTA,
        FASTA_INDEX = rules.qc_samtools_index.output.INDEX_GENOME_PREFIX
    output:
        VCF =  Path(config['working_dir']) / 'qc' / '{name}' / 'freebayes' / 'variants.vcf',
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / 'freebayes' / 'commands.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'freebayes',
        freebayes_options = config.get('freebayes',{})
    run:
        from camel.app.tools.freebayes.freebayes import Freebayes
        dir_working = Path(str(params.running_dir)).absolute()
        freebayes = Freebayes(camel)
        freebayes.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))], 'FASTA': [ToolIOFile(Path(input.FASTA))]})
        freebayes.update_parameters(**params.freebayes_options)
        step = Step(str(rule), freebayes, camel, dir_working, config)
        step.run_step()
        freebayes.informs['_tag'] = wildcards.name
        SnakemakeUtils.dump_object(freebayes.informs, Path(output.INFORMS))

rule qc_sniffles:
    """
    Checks for structural variants in the assembly using long reads.
    """
    input:
        BAM = rules.qc_bam_sorting_ont.output.BAM,
        BAM_INDEX = rules.qc_bam_indexing_ont.output.BAM,
        FASTA = rules.qc_copy_fasta_file.output.FASTA,
        FASTA_INDEX = rules.qc_samtools_index.output.INDEX_GENOME_PREFIX
    output:
        VCF = Path(config['working_dir']) / 'qc' / '{name}' / 'sniffles' / 'variants.vcf',
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / 'sniffles' / 'commands.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'sniffles'
    threads: 4
    run:
        from camel.app.tools.sniffles.sniffles import Sniffles
        dir_working = Path(str(params.running_dir)).absolute()
        sniffles = Sniffles(camel)
        sniffles.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))], 'FASTA': [ToolIOFile(Path(input.FASTA))]})
        sniffles.update_parameters(threads=threads)
        step = Step(str(rule), sniffles, camel, dir_working, config)
        step.run_step()
        sniffles.informs['_tag'] = wildcards.name
        SnakemakeUtils.dump_object(sniffles.informs, Path(output.INFORMS))

rule qc_clair3:
    """
    Checks for small variants in the assembly using long reads.
    """
    input:
        BAM = rules.qc_bam_sorting_ont.output.BAM,
        BAM_INDEX = rules.qc_bam_indexing_ont.output.BAM,
        FASTA = rules.qc_copy_fasta_file.output.FASTA,
        FASTA_INDEX = rules.qc_samtools_index.output.INDEX_GENOME_PREFIX
    output:
        VCF_GZ = Path(config['working_dir']) / 'qc' / '{name}' / 'clair3_output' / 'merge_output.vcf.gz',
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / 'clair3_output' / 'commands.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}',
        clair3_options = config.get('clair3', {})
    threads: 8
    run:
        from camel.app.tools.clair3.clair3 import Clair3
        dir_working = Path(str(params.running_dir)).absolute()
        clair3 = Clair3(camel)
        clair3.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))], 'FASTA': [ToolIOFile(Path(input.FASTA))]})
        clair3.update_parameters(
            **params.clair3_options, chunk_size=100_000, platform='ont', no_phasing=True, include_ctgs=True,
            threads=threads)
        step = Step(str(rule), clair3, camel, dir_working, config)
        step.run_step()
        clair3.informs['_tag'] = wildcards.name
        SnakemakeUtils.dump_object(clair3.informs, Path(output.INFORMS))

rule qc_unzip_clair3_vcf:
    """
    Unzips the output from clair3 to be compatible with the qc_add_vcf_info_to_informs rule.
    """
    input:
        VCF_GZ = rules.qc_clair3.output.VCF_GZ
    output:
        VCF = Path(config['working_dir']) / 'qc' / '{name}' / 'clair3_output' / 'variants.vcf'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'clair3_output'
    run:
        from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
        dir_working = Path(str(params.running_dir)).absolute()
        bcftools_view = BcftoolsView(camel)
        bcftools_view.add_input_files({'VCF_GZ': [ToolIOFile(Path(input.VCF_GZ))]})
        bcftools_view.update_parameters(output_format='VCF', compress_output=False, output_filename='variants.vcf')
        step = Step(str(rule), bcftools_view, camel, dir_working, config)
        step.run_step()

rule qc_add_vcf_info_to_informs:
    """
    Parses the VCF file and adds the statistics to the informs.
    """
    input:
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / '{method}' / 'commands.io',
        VCF = Path(config['working_dir']) / 'qc' / '{name}' / '{method}' / 'variants.vcf'
    output:
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / '{method}' / 'informs.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / wildcards.name / wildcards.method
    run:
        from camel.app.components.vcf.vcfutils import VCFUtils
        informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        if wildcards.method == 'sniffles':
            informs['nb_of_variants'] = sum(informs['variants'].values())
            informs['nb_of_snps'] = 0
            informs['nb_of_indels'] = informs['variants']['DEL'] + informs['variants']['INS']
            informs['nb_of_svs'] = informs['variants']['BND'] + informs['variants']['DUP'] + informs['variants']['INV']
        else:
            variants = VCFUtils.retrieve_variants(Path(input.VCF))
            informs['nb_of_variants'] = len(variants)
            informs['nb_of_snps'] = sum(v.is_snp for v in variants)
            informs['nb_of_indels'] = sum(v.is_indel for v in variants)
            informs['nb_of_svs'] = 0
        SnakemakeUtils.dump_object(informs, Path(output.INFORMS))

rule qc_ale:
    """
    Generates ALE QC report for the final assembly.
    """
    input:
        SAM = rules.qc_read_mapping_illumina.output.SAM,
        FASTA = rules.qc_copy_fasta_file.output.FASTA,
        FASTA_INDEX = rules.qc_samtools_index.output.INDEX_GENOME_PREFIX
    output:
        ALE = Path(config['working_dir']) / 'qc' / '{name}' / 'ale_illumina' / 'ALE.ale',
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / 'ale_illumina' / 'informs-report.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'ale_illumina'
    run:
        from camel.app.tools.ale.ale import ALE
        dir_working = Path(str(params.running_dir)).absolute()
        ale_report = ALE(camel)
        ale_report.add_input_files({'SAM': [ToolIOFile(Path(input.SAM))], 'FASTA': [ToolIOFile(Path(input.FASTA))]})
        step = Step(str(rule), ale_report, camel, dir_working, config)
        step.run_step()
        ale_report.informs['_tag'] = f'{wildcards.name}'
        SnakemakeUtils.dump_object(ale_report.informs, Path(output.INFORMS))

rule qc_ale2wiggle:
    """
    Generates wiggle files from ALE for IGV visualization.
    """
    input:
        ALE = rules.qc_ale.output.ALE
    output:
        TSV_1 = Path(config['working_dir']) / 'qc' / '{name}' / 'ale_illumina' / 'ALE.ale-depth.wig',
        TSV_2 = Path(config['working_dir']) / 'qc' / '{name}' / 'ale_illumina' / 'ALE.ale-kmer.wig',
        TSV_3 = Path(config['working_dir']) / 'qc' / '{name}' / 'ale_illumina' / 'ALE.ale-insert.wig',
        TSV_4 = Path(config['working_dir']) / 'qc' / '{name}' / 'ale_illumina' / 'ALE.ale-place.wig',
        INFORMS = Path(config['working_dir']) / 'qc' / '{name}' / 'ale_illumina' / 'commands-wiggle.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'qc' / f'{wildcards.name}' / 'ale_illumina'
    run:
        from camel.app.tools.ale.ale2wiggle import ALE2Wiggle
        dir_working = Path(str(params.running_dir)).absolute()
        ale2wiggle_report = ALE2Wiggle(camel)
        ale2wiggle_report.add_input_files({'ALE': [ToolIOFile(Path(input.ALE))]})
        step = Step(str(rule), ale2wiggle_report, camel, dir_working, config)
        step.run_step()
        ale2wiggle_report.informs['_tag'] = f'{wildcards.name}'
        SnakemakeUtils.dump_object(ale2wiggle_report.informs, Path(output.INFORMS))
