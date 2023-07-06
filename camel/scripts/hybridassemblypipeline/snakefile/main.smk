import gzip
import shutil
from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import assembly_flye, medaka_polishing, short_read_polishing
from camel.scripts.hybridassemblypipeline.snakefile import quality_checks

camel = Camel.get_instance()

include: assembly_flye.SNAKEFILE_ASSEMBLY_FLYE
include: medaka_polishing.SNAKEFILE_MEDAKA_POLISHING
include: short_read_polishing.SNAKEFILE_POLISHING
include: quality_checks.SNAKEFILE_QC

assembly_steps = ['Flye', 'Medaka', 'Polypolish', 'POLCA', 'Unicycler']

#########
# Rules #
#########
rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        Path(config['working_dir']) / config['output_html']

rule trim_illumina:
    """
    Trims the Illumina reads.
    """
    input:
        FQ_fwd = config['input']['illumina'][0],
        FQ_rev = config['input']['illumina'][1]
    output:
        FQ_dict = Path(config['working_dir']) / 'trimming' / 'illumina' / 'fq_dict.io',
        TSV = Path(config['working_dir']) / 'trimming' / 'illumina' / 'trimming_illumina.tsv',
        HTML = Path(config['working_dir']) / 'trimming' / 'illumina' / 'html.io'
    params:
        dir_ = Path(config['working_dir']) / 'trimming' / 'illumina'
    threads: 4
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.components.workflows.trimmingilluminawrapper import TrimmingIlluminaWrapper
        wrapper = TrimmingIlluminaWrapper(Path(params.dir_))
        wrapper.run_workflow([Path(input.FQ_fwd), Path(input.FQ_rev)], threads=threads)
        wrapper.output.tsv_summary.rename(Path(output.TSV))
        SnakemakeUtils.dump_object(wrapper.output.report_section, Path(output.HTML))
        fq_out = FastqInput(
            read_type='illumina',
            pe=wrapper.output.trimmed_reads_pe,
            se_fwd=wrapper.output.trimmed_reads_se_fwd,
            se_rev=wrapper.output.trimmed_reads_se_rev,
            is_trimmed=True,
            is_pe=True)
        SnakemakeUtils.dump_object(fq_out.to_fq_dict(), Path(output.FQ_dict))

rule trim_ont_workflow:
    """
    Trims the ONT reads.
    """
    input:
        FASTQ = config['input']['ont']
    output:
        FASTQ = Path(config['working_dir']) / 'trimming' / 'ont' / 'fastq.io',
        HTML = Path(config['working_dir']) / 'trimming' / 'ont' / 'html.io'
    params:
        dir_ = Path(config['working_dir']) / 'trimming' / 'ont'
    threads: 4
    run:
        from camel.app.components.workflows.trimmingontwrapper import TrimmingONTWrapper
        wrapper = TrimmingONTWrapper(Path(params.dir_).absolute())
        wrapper.run_workflow(Path(input.FASTQ), threads=threads)
        wrapper.output.trimmed_reads[0].path.rename(Path(output.FASTQ))
        SnakemakeUtils.dump_object(wrapper.output.report_section, Path(output.HTML))

rule set_trimming_ont_output:
    """
    This rule gzip the filtlong output reads into the correct location.
    """
    input:
        FASTQ = rules.trim_ont_workflow.output.FASTQ
    output:
        FASTQ = Path(config['working_dir']) / 'trimming' / 'ont' / 'trimmed.fastq.gz',
        FASTQ_IO = Path(config['working_dir']) / 'fq_dict.io'
    params:
        dir_ = Path(config['working_dir']) / 'trimming' / 'ont'
    run:
        input_fastq = open(input.FASTQ, 'rb').read()
        with gzip.open(output.FASTQ, 'wb') as handle:
            handle.write(input_fastq)
        SnakemakeUtils.dump_object({'SE': [ToolIOFile(Path(output.FASTQ))]}, Path(output.FASTQ_IO))

rule unicycler:
    """
    Runs unicycler, which is a short-read first approach to assemble reads.
    """
    input:
        FQ_dict = rules.trim_illumina.output.FQ_dict,
        FQ_SE = rules.trim_ont_workflow.output.FASTQ
    output:
        FASTA = Path(config['working_dir']) / 'unicycler' / 'assembly.fasta',
        INFORMS = Path(config['working_dir']) / 'unicycler' / 'commands.io'
    params:
        working_dir = Path(config['working_dir'])
    threads: 16
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.unicycler.unicycler import Unicycler
        unicycler_assembly = Unicycler(camel)
        fq_illumina = FastqInput.from_fq_dict(Path(input.FQ_dict), 'illumina')
        unicycler_assembly.add_input_files({
            'FASTQ_PE': fq_illumina.pe,
            'FASTQ_SE': [ToolIOFile(Path(input.FQ_SE))]})
        unicycler_assembly.update_parameters(output_dir='unicycler', threads=threads)
        step = Step(str(rule), unicycler_assembly, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_object(unicycler_assembly.informs, Path(output.INFORMS))

rule copy_medaka_to_short_read_polishing:
    input:
        FASTA = Path(config['working_dir']) / medaka_polishing.OUTPUT_ASSEMBLY_FASTA,
        FASTQ = rules.trim_illumina.output.FQ_dict
    output:
        FASTA = Path(config['working_dir']) / short_read_polishing.INPUT_ASSEMBLY_FASTA,
        FASTQ = Path(config['working_dir']) / short_read_polishing.INPUT_READS_FASTQ
    run:
        shutil.copyfile(input.FASTA, output.FASTA)
        shutil.copyfile(input.FASTQ, output.FASTQ)

rule combine_informs_quast:
    """
    Combines the QUAST informs into a single table
    """
    input:
        INFORMS_quast = [Path(config['working_dir']) / 'qc' / name / 'quast' / 'informs.io' for name in assembly_steps]
    output:
        TSV = Path(config['working_dir'] / 'report' / 'quast_stats-all.tsv')
    run:
        records_out = []
        for path_informs in [Path(x) for x in input.INFORMS_quast]:
            quast_informs = SnakemakeUtils.load_object(path_informs)
            assembly_key = path_informs.parent.parent.name
            records_out.append({
                'Assembly step': assembly_key,
                'N50': '{:,}'.format(int(quast_informs['contig']['N50'])),
                'Nb. of contigs': '{:,}'.format(int(quast_informs['contig']['# contigs (>= 1000 bp)'])),
                'Total length': '{:,}'.format(int(quast_informs['genome']['Total length']))
            })
        pd.DataFrame(records_out).to_csv(output.TSV, sep='\t', index=False)

rule combine_informs_variant_calling_short_reads:
    """
    Combines the informs from the variant calling with short reads.
    """
    input:
        INFORMS_freebayes = [Path(config['working_dir']) / 'qc' / name / 'freebayes' / 'informs.io' for name in assembly_steps],
        INFORMS_clair3 = [Path(config['working_dir']) / 'qc' / name / 'clair3_output' / 'informs.io' for name in assembly_steps]
    output:
        TSV = Path(config['working_dir'] / 'report' / 'variant_calling_all-short.tsv')
    run:
        records_out = []
        for path_freebayes, path_clair3 in zip(
                [Path(x) for x in input.INFORMS_freebayes], [Path(x) for x in input.INFORMS_clair3]):
            informs_freebayes = SnakemakeUtils.load_object(path_freebayes)
            informs_clair3 = SnakemakeUtils.load_object(path_clair3)
            assembly_key = path_freebayes.parent.parent.name
            records_out.append({
                'Assembly step': assembly_key,
                'Nb. of SNPs (FreeBayes)': '{:,}'.format(int(informs_freebayes['nb_of_snps'])),
                'Nb. of indels (FreeBayes)': '{:,}'.format(int(informs_freebayes['nb_of_indels'])),
                'Nb. of SNPs (Clair3)': '{:,}'.format(int(informs_clair3['nb_of_snps'])),
                'Nb. of indels (Clair3)': '{:,}'.format(int(informs_clair3['nb_of_indels']))
            })
        pd.DataFrame(records_out).to_csv(output.TSV, sep='\t', index=False)

rule report_long_variant_calling:
    """
    Rule to generate the report table for sniffles.
    """
    input:
        INFORMS_sniffles = [Path(config['working_dir']) / 'qc' / name / 'sniffles' / 'informs.io' for name in assembly_steps]
    output:
        TSV = Path(config['working_dir'] / 'report' / 'variant_calling-sniffles.tsv')
    run:
        sniffles_table = []
        for path_informs in [Path(x) for x in input.INFORMS_sniffles]:
            informs_sniffles = SnakemakeUtils.load_object(path_informs)
            assembly_key = path_informs.parent.parent.name
            sniffles_table.append({
                'Assembly step': assembly_key,
                'Nb. of indels': '{:,}'.format(int(informs_sniffles['nb_of_indels'])),
                'Nb. of SVs': '{:,}'.format(int(informs_sniffles['nb_of_svs'])),
            })
        pd.DataFrame(sniffles_table).to_csv(output.TSV, sep='\t', index=False)

rule report_mapping_stats:
    """
    Rule to generate the report table for the mapping statistics.
    """
    input:
        INFORMS_mapping_illumina = [Path(config['working_dir']) / 'qc' / name / 'read_mapping' / 'illumina' / 'flagstat.io' for name in assembly_steps],
        INFORMS_mapping_ont = [Path(config['working_dir']) / 'qc' / name / 'read_mapping' / 'ont' / 'flagstat-longreads.io'  for name in assembly_steps],
        INFORMS_depth_illumina = [Path(config['working_dir']) / 'qc' / name / 'read_mapping' / 'illumina' / 'samtools-depth.io' for name in assembly_steps],
        INFORMS_depth_ont = [Path(config['working_dir']) / 'qc' / name / 'read_mapping' / 'ont' / 'samtools-depth-long.io' for name in assembly_steps]
    output:
        TSV = Path(config['working_dir'] / 'report' / 'mapping_statistics.tsv')
    run:
        records_out = []
        for path_map_illumina, path_map_ont, path_depth_illumina, path_depth_ont in zip(
                [Path(x) for x in input.INFORMS_mapping_illumina], [Path(x) for x in input.INFORMS_mapping_ont],
                [Path(x) for x in input.INFORMS_depth_illumina], [Path(x) for x in input.INFORMS_depth_ont]):
            assembly_key = path_map_illumina.parent.parent.parent.name
            illumina_mapping_informs = SnakemakeUtils.load_object(path_map_illumina)
            ont_mapping_informs = SnakemakeUtils.load_object(path_map_ont)
            illumina_depth_informs = SnakemakeUtils.load_object(path_depth_illumina)
            ont_depth_informs = SnakemakeUtils.load_object(path_depth_ont)
            records_out.append({
                'Assembly step': assembly_key,
                'Mapping rate (Illumina)': '{:.2f}'.format(illumina_mapping_informs['mapped'][0] / illumina_mapping_informs['total'][0] * 100),
                'Mapping rate (ONT)': '{:.2f}'.format(ont_mapping_informs['mapped'][0] / ont_mapping_informs['total'][0] * 100),
                'Median depth (Illumina)': int(illumina_depth_informs['median_depth']),
                'Median depth (ONT)': int(ont_depth_informs['median_depth'])
            })
        pd.DataFrame(records_out).to_csv(output.TSV,sep='\t',index=False)

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
    Creates a report section with the commands used in the pipeline. 
    """
    input:
        unicycler_commands = Path(config['working_dir']) / 'unicycler' / 'commands.io',
        flye_commands = Path(config['working_dir']) / 'assembly_flye' / 'flye' / 'informs.io',
        medaka_consensus_commands = Path(config['working_dir']) / 'medaka' / 'commands-consensus.io',
        medaka_stitch_commands = Path(config['working_dir']) / 'medaka' / 'commands-stitch.io',
        polypolish_commands = Path(config['working_dir']) / 'polishing' / 'polypolish'  / 'polypolish.io',
        polca_commands = Path(config['working_dir']) / 'polishing' / 'polca' / 'polca.io',
        quast_commands = [Path(config['working_dir']) / 'qc' / name / 'quast' / 'commands.io' for name in assembly_steps],
        quast_combined_commands = Path(config['working_dir']) / 'qc' / 'quast_combined' / 'commands.io',
        bwa_commands = [Path(config['working_dir']) / 'qc' / name / 'read_mapping' / 'illumina' / 'commands.io' for name in assembly_steps],
        freebayes_commands = [Path(config['working_dir']) / 'qc' / name / 'freebayes' / 'commands.io' for name in assembly_steps],
        sniffles_commands = [Path(config['working_dir']) / 'qc' / name / 'sniffles' / 'commands.io' for name in assembly_steps],
        clair3_commands = [Path(config['working_dir']) / 'qc' / name / 'clair3_output' / 'commands.io' for name in assembly_steps],
        ale_report_commands = [Path(config['working_dir']) / 'qc' / name / 'ale_illumina' / 'informs-report.io' for name in assembly_steps],
        ale_wiggle_commands = [Path(config['working_dir']) / 'qc' / name / 'ale_illumina' / 'commands-wiggle.io' for name in assembly_steps]
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        working_dir = config['working_dir']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        informs = []
        for content in [SnakemakeUtils.load_object(Path(io)) for io in input]:
            if isinstance(content, dict):
                informs.append(content)
            elif isinstance(content, list):
                informs.extend(content)
        section = SnakePipelineUtils.create_commands_section(informs, params.working_dir)
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule report_create_sections:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        FASTA = [Path(config['working_dir']) / 'qc' / key / 'consensus.fasta' for key in quality_checks.consensus_by_tool.keys()],
        TSV_quast = rules.combine_informs_quast.output.TSV,
        TSV_vc = rules.combine_informs_variant_calling_short_reads.output.TSV,
        TSV_mapping = rules.report_mapping_stats.output.TSV,
        TSV_sniffles = rules.report_long_variant_calling.output.TSV,
        report_citations = rules.report_pickle_citations.output.HTML,
        INFORMS_quast = Path(config['working_dir']) / 'qc' / assembly_steps[0] / 'quast' / 'commands.io',
        INFORMS_freebayes = Path(config['working_dir']) / 'qc' / assembly_steps[0] / 'freebayes' / 'commands.io',
        INFORMS_clair3 = Path(config['working_dir']) / 'qc' / assembly_steps[0] / 'clair3_output' / 'commands.io',
        INFORMS_sniffles = [Path(config['working_dir']) / 'qc' / steps / 'sniffles' / 'informs.io' for steps in assembly_steps],
        INFORMS_mapping = Path(config['working_dir']) / 'qc' / assembly_steps[0] / 'read_mapping' / 'illumina' / 'commands.io',
        INFORMS_ale = [Path(config['working_dir']) / 'qc' / steps / 'ale_illumina' / 'informs-report.io' for steps in assembly_steps],
        report_commands = rules.report_command_section.output.HTML,
        HTML_trim_illumina = rules.trim_illumina.output.HTML,
        HTML_trim_ont = rules.trim_ont_workflow.output.HTML,
        HTML_quast = Path(config['working_dir']) / 'qc' / 'quast_combined' / 'report.html',
        VCF_sniffles = [Path(config['working_dir']) / 'qc' / key / 'sniffles' / 'variants.vcf' for key in quality_checks.consensus_by_tool.keys()],
        WIGGLE_ale = [Path(config['working_dir']) / 'qc' / key / 'ale_illumina' / f'ALE.ale-{ale_key}.wig' for key in quality_checks.consensus_by_tool.keys() for ale_key in quality_checks.ALE_KEYS]
    output:
        HTML = Path(config['working_dir']) / Path(config['output_html'])
    params:
        sample_name = config['sample_name'],
        input = config['input'],
        pipeline = config['pipeline'],
        working_dir = Path(config['working_dir']),
        output_dir = config['output_dir']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.scripts.hybridassemblypipeline.reporter.hybridassemblyreporter import HybridAssemblyReporter
        reporter = HybridAssemblyReporter(camel)
        reporter.add_input_files({
            'FASTA': [ToolIOFile(Path(p)) for p in input.FASTA],
            'TSV_quast': [ToolIOFile(Path(input.TSV_quast))],
            'TSV_vc': [ToolIOFile(Path(input.TSV_vc))],
            'TSV_sniffles': [ToolIOFile(Path(input.TSV_sniffles))],
            'VCF_sniffles': [ToolIOFile(Path(x)) for x in input.VCF_sniffles],
            'TSV_mapping': [ToolIOFile(Path(input.TSV_mapping))],
            'HTML_trim_illumina': [ToolIOValue(SnakemakeUtils.load_object(Path(input.HTML_trim_illumina)))],
            'HTML_trim_ont': [ToolIOValue(SnakemakeUtils.load_object(Path(input.HTML_trim_ont)))],
            'HTML_quast': [ToolIOFile(Path(input.HTML_quast))],
            'WIGGLE_ale': [ToolIOFile(Path(x)) for x in input.WIGGLE_ale]
        })
        reporter.add_input_informs({
            'quast': SnakemakeUtils.load_object(Path(input.INFORMS_quast)),
            'freebayes': SnakemakeUtils.load_object(Path(input.INFORMS_freebayes)),
            'clair3': SnakemakeUtils.load_object(Path(input.INFORMS_clair3)),
            'sniffles': [SnakemakeUtils.load_object(Path(f)) for f in input.INFORMS_sniffles],
            'mapping': SnakemakeUtils.load_object(Path(input.INFORMS_mapping)),
            'commands': SnakemakeUtils.load_object(Path(input.report_commands)),
            'citations': SnakemakeUtils.load_object(Path(input.report_citations)),
            'ale': [SnakemakeUtils.load_object(Path(f)) for f in input.INFORMS_ale],
            'sample_name': params.sample_name,
            'pipeline': params.pipeline,
            'input': params.input
        })
        reporter.update_parameters(output_filename=str(output.HTML), output_dir=str(params.output_dir))
        step = Step(str(rule), reporter, camel, params.working_dir, config)
        step.run_step()
