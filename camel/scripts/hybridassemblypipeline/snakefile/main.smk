import itertools
import shutil
from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import core, downsampling, polish_assembly_long, polish_assembly_short, \
    trimming_illumina, trimming_ont, assembly
from camel.scripts.hybridassemblypipeline.snakefile import qc_hybrid

camel = Camel.get_instance()


#######################
# Included Snakefiles #
#######################
include: core.SNAKEFILE_CORE
include: downsampling.SNAKEFILE_DOWNSAMPLING
include: trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA
include: trimming_ont.SNAKEFILE_TRIMMING_ONT
include: assembly.SNAKEFILE_ASSEMBLY
include: qc_hybrid.SNAKEFILE_QC

ruleorder: copy_medaka_unicycler_to_short_read_polishing > core_link_fasta_to_polishing

#########
# Rules #
#########
rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        Path(config['working_dir']) / config['output_html']

rule unicycler:
    """
    Runs unicycler, which is a short-read first approach to assemble reads.
    """
    input:
        FQ_dict = Path(config['working_dir']) / 'fq_dict.io'
    output:
        FASTA = Path(config['working_dir']) / 'unicycler' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / 'unicycler' / 'commands.io'
    params:
        working_dir = Path(config['working_dir'])
    threads: 16
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.unicycler.unicycler import Unicycler
        unicycler_assembly = Unicycler(camel)
        fq_hybrid = FastqInput.from_fq_dict(Path(input.FQ_dict), 'hybrid')
        unicycler_assembly.add_input_files({
            'FASTQ_PE': fq_hybrid.pe,
            'FASTQ_SE': fq_hybrid.se})
        unicycler_assembly.update_parameters(output_dir='unicycler', threads=threads)
        step = Step(str(rule), unicycler_assembly, camel, params.working_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(unicycler_assembly, output)

rule copy_unicycler_assembly_to_medaka_input:
    input:
        FASTA_unicycler = Path(config['working_dir']) / 'unicycler' / 'fasta.io'
    output:
        FASTA_medaka_unicycler = str(Path(config['working_dir']) / polish_assembly_long.INPUT_ASSEMBLY_FASTA).format(assembly_type='unicycler')
    run:
        shutil.copyfile(input.FASTA_unicycler, output.FASTA_medaka_unicycler)

rule copy_medaka_unicycler_to_short_read_polishing:
    input:
        FASTA_unicycler = str(Path(config['working_dir']) / polish_assembly_long.OUTPUT_ASSEMBLY_FASTA).format(assembly_type='unicycler')
    output:
        FASTA_unicycler = str(Path(config['working_dir']) / polish_assembly_short.INPUT_ASSEMBLY_FASTA).format(assembly_type='unicycler')
    run:
        shutil.copyfile(input.FASTA_unicycler, output.FASTA_unicycler)

rule combine_informs_quast:
    """
    Combines the QUAST informs into a single table
    """
    input:
        INFORMS_quast = [Path(config['working_dir']) / 'qc_hybrid' / name / 'quast' / 'informs.io' for name in config['assembly_steps']]
    output:
        TSV = Path(config['working_dir']) / 'report' / 'quast_stats-all-tsv.tsv'
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
        INFORMS_freebayes = [Path(config['working_dir']) / 'qc_hybrid' / name / 'freebayes' / 'informs.io' for name in config['assembly_steps']],
        INFORMS_clair3 = [Path(config['working_dir']) / 'qc_hybrid' / name / 'clair3_output' / 'informs.io' for name in config['assembly_steps']]
    output:
        TSV = Path(config['working_dir']) / 'report' / 'variant_calling_all-short.tsv'
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
        INFORMS_sniffles = [Path(config['working_dir']) / 'qc_hybrid' / name / 'sniffles' / 'informs.io' for name in config['assembly_steps']]
    output:
        TSV = Path(config['working_dir']) / 'report' / 'variant_calling-sniffles.tsv'
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
        INFORMS_mapping_illumina = [Path(config['working_dir']) / 'qc_hybrid' / name / 'read_mapping' / 'illumina' / 'flagstat.io' for name in config['assembly_steps']],
        INFORMS_mapping_ont = [Path(config['working_dir']) / 'qc_hybrid' / name / 'read_mapping' / 'ont' / 'flagstat-longreads.io'  for name in config['assembly_steps']],
        INFORMS_depth_illumina = [Path(config['working_dir']) / 'qc_hybrid' / name / 'read_mapping' / 'illumina' / 'samtools-depth.io' for name in config['assembly_steps']],
        INFORMS_depth_ont = [Path(config['working_dir']) / 'qc_hybrid' / name / 'read_mapping' / 'ont' / 'samtools-depth-long.io' for name in config['assembly_steps']]
    output:
        TSV = Path(config['working_dir']) / 'report' / 'mapping_statistics.tsv'
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
        pd.DataFrame(records_out).to_csv(output.TSV, sep='\t', index=False)

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
        unicycler_commands = Path(config['working_dir']) / 'unicycler' / 'commands.io' if 'unicycler' in config['base_assemblies'] else [],
        flye_commands = assembly.get_command_informs(config),
        medaka_inference_commands = [Path(config['working_dir']) / 'polish' / 'long_reads' / name / 'inference' / 'commands-inference.io' for name in config['base_assemblies']],
        medaka_sequence_commands = [Path(config['working_dir']) / 'polish' / 'long_reads' / name / 'sequence' / 'commands-sequence.io' for name in config['base_assemblies']],
        polypolish_commands = [Path(config['working_dir']) / 'polish' / 'short_reads' / name / 'polypolish' / 'informs.io' for name in config['base_assemblies']],
        polca_commands = [Path(config['working_dir']) / 'polish' / 'short_reads' / name / 'polca' / 'informs.io' for name in config['base_assemblies']],
        quast_commands = [Path(config['working_dir']) / 'qc_hybrid' / name / 'quast' / 'commands.io' for name in config['assembly_steps']],
        quast_combined_commands = Path(config['working_dir']) / 'qc_hybrid' / 'quast_combined' / 'commands.io',
        bwa_commands = [Path(config['working_dir']) / 'qc_hybrid' / name / 'read_mapping' / 'illumina' / 'commands.io' for name in config['assembly_steps']],
        freebayes_commands = [Path(config['working_dir']) / 'qc_hybrid' / name / 'freebayes' / 'commands.io' for name in config['assembly_steps']],
        sniffles_commands = [Path(config['working_dir']) / 'qc_hybrid' / name / 'sniffles' / 'commands.io' for name in config['assembly_steps']],
        clair3_commands = [Path(config['working_dir']) / 'qc_hybrid' / name / 'clair3_output' / 'commands.io' for name in config['assembly_steps']],
        ale_report_commands = [Path(config['working_dir']) / 'qc_hybrid' / name / 'ale_illumina' / 'informs-report.io' for name in config['assembly_steps']],
        ale_wiggle_commands = [Path(config['working_dir']) / 'qc_hybrid' / name / 'ale_illumina' / 'commands-wiggle.io' for name in config['assembly_steps']]
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.export_command_section(input, Path(output.HTML), Path(params.dir_))

rule report_create_sections:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        FASTA = [Path(config['working_dir']) / 'qc_hybrid' / key / 'fasta-index.io' for key in config['assembly_steps']],
        TSV_quast = rules.combine_informs_quast.output.TSV,
        TSV_vc = rules.combine_informs_variant_calling_short_reads.output.TSV,
        TSV_mapping = rules.report_mapping_stats.output.TSV,
        TSV_sniffles = rules.report_long_variant_calling.output.TSV,
        report_citations = rules.report_pickle_citations.output.HTML,
        INFORMS_quast = Path(config['working_dir']) / 'qc_hybrid' / config['assembly_steps'][0] / 'quast' / 'commands.io',
        INFORMS_freebayes = Path(config['working_dir']) / 'qc_hybrid' / config['assembly_steps'][0] / 'freebayes' / 'commands.io',
        INFORMS_clair3 = Path(config['working_dir']) / 'qc_hybrid' / config['assembly_steps'][0] / 'clair3_output' / 'commands.io',
        INFORMS_sniffles = [Path(config['working_dir']) / 'qc_hybrid' / steps / 'sniffles' / 'informs.io' for steps in config['assembly_steps']],
        INFORMS_mapping = Path(config['working_dir']) / 'qc_hybrid' / config['assembly_steps'][0] / 'read_mapping' / 'illumina' / 'commands.io',
        INFORMS_ale = [Path(config['working_dir']) / 'qc_hybrid' / steps / 'ale_illumina' / 'informs-report.io' for steps in config['assembly_steps']],
        report_commands = rules.report_command_section.output.HTML,
        HTML_trimming_illumina = Path(config['working_dir']) / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_REPORT,
        HTML_trimming_ont = Path(config['working_dir']) / trimming_ont.OUTPUT_TRIMMING_ONT_REPORT,
        HTML_quast = Path(config['working_dir']) / 'qc_hybrid' / 'quast_combined' / 'report.html',
        VCF_sniffles = [Path(config['working_dir']) / 'qc_hybrid' / key / 'sniffles' / 'vcf.io' for key in config['assembly_steps']],
        WIGGLE_ale = [Path(config['working_dir']) / 'qc_hybrid' / key / 'ale_illumina' / 'wiggle.io' for key in config['assembly_steps']]
    output:
        HTML = Path(config['working_dir']) / Path(config['output_html'])
    params:
        sample_name = config['sample_name'],
        input = config['input'],
        pipeline = config['pipeline'],
        working_dir = Path(config['working_dir']),
        output_dir = config['output_dir']
    run:
        from camel.scripts.hybridassemblypipeline.reporter.hybridassemblyreporter import HybridAssemblyReporter
        reporter = HybridAssemblyReporter(camel)
        reporter.add_input_files({
            'FASTA': [SnakemakeUtils.load_object(Path(FASTA))[0] for FASTA in input.FASTA],
            'TSV_quast': [ToolIOFile(Path(input.TSV_quast))],
            'TSV_vc': [ToolIOFile(Path(input.TSV_vc))],
            'TSV_sniffles': [ToolIOFile(Path(input.TSV_sniffles))],
            'VCF_sniffles': [SnakemakeUtils.load_object(Path(VCF))[0] for VCF in input.VCF_sniffles],
            'TSV_mapping': [ToolIOFile(Path(input.TSV_mapping))],
            'HTML_trim_illumina': [SnakemakeUtils.load_object(Path(input.HTML_trimming_illumina))[0]],
            'HTML_trim_ont': [SnakemakeUtils.load_object(Path(input.HTML_trimming_ont))[0]],
            'HTML_quast': [ToolIOFile(Path(input.HTML_quast))],
            'WIGGLE_ale': list(itertools.chain(*[SnakemakeUtils.load_object(Path(x))[:] for x in input.WIGGLE_ale]))
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
        step = Step(str(rule), reporter, camel, params.working_dir)
        step.run_step()
