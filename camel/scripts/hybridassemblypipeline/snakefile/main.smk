import gzip
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.hybridassemblypipeline.snakefile import assembly_flye, short_read_polishing, medaka_snakemake, quality_checks

camel = Camel.get_instance()

include: assembly_flye.SNAKEFILE_FLYE
include: medaka_snakemake.SNAKEFILE_POLISHING
include: short_read_polishing.SNAKEFILE_POLISHING
include: quality_checks.SNAKEFILE_QC

assembly_steps = ['Flye', 'Medaka', 'POLCA', 'Polypolish', 'Unicycler']

#########
# Rules #
#########

rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        Path(config['working_dir']) / config['output_html'],
        Path(config['working_dir'] / config['output'])

rule trim_illumina:
    """
    This rule trims the illumina reads using trimmomatic.
    """
    input:
        FQ_fwd = config['input']['illumina'][0],
        FQ_rev = config['input']['illumina'][1]
    output:
        FQ_1P = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_1P.fastq.gz",
        FQ_2P = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_2P.fastq.gz",
        FQ_1S = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_1U.fastq.gz",
        FQ_2S = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_2U.fastq.gz",
        TSV = Path(config['working_dir']) / 'trimming' / 'illumina' / 'trimming_illumina.tsv',
        HTML = Path(config['working_dir']) / 'trimming' / 'illumina' / 'html.io'
    params:
        dir_ = Path(config['working_dir']) / 'trimming' / 'illumina'
    threads: 4
    run:
        from camel.app.components.workflows.trimmingilluminawrapper import TrimmingIlluminaWrapper
        wrapper = TrimmingIlluminaWrapper(Path(params.dir_).absolute())
        wrapper.run_workflow([Path(input.FQ_fwd), Path(input.FQ_rev)], threads=threads)
        wrapper.output.trimmed_reads_pe[0].path.rename(Path(output.FQ_1P))
        wrapper.output.trimmed_reads_pe[1].path.rename(Path(output.FQ_2P))
        wrapper.output.trimmed_reads_se_fwd[0].path.rename(Path(output.FQ_1S))
        wrapper.output.trimmed_reads_se_rev[0].path.rename(Path(output.FQ_2S))
        wrapper.output.tsv_summary.rename(Path(output.TSV))
        SnakemakeUtils.dump_object(wrapper.output.report_section, Path(output.HTML))

rule trim_ont_workflow:
    """
    This rule trims the ONT reads using the trimmingont wrapper.
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
        wrapper.run_workflow(Path(input.FASTQ), threads=threads, )
        wrapper.output.trimmed_reads[0].path.rename(Path(output.FASTQ))
        SnakemakeUtils.dump_object(wrapper.output.report_section, Path(output.HTML))

# rule trim_ont:
#     """
#     This rule trims the ONT reads using filtlong.
#     """
#     input:
#         FASTQ = config['input']['ont']
#     output:
#         FASTQ = Path(config['working_dir']) / 'trimming' / 'ont' / 'fastq.io'
#     params:
#         working_directory = Path(config['working_dir']) / 'trimming' / 'ont',
#         filtlong_options= config.get('filtlong',{})
#     threads: 4
#     run:
#         from camel.app.tools.filtlong.filtlong import Filtlong
#         filtlong = Filtlong(camel)
#         filtlong.add_input_files({'FASTQ': [ToolIOFile(Path(input.FASTQ))]})
#         filtlong.update_parameters(**params.filtlong_options)
#         step = Step(str(rule), filtlong, camel, params.working_directory, config)
#         step.run_step()
#         SnakemakeUtils.dump_tool_outputs(filtlong, output)

rule set_trimming_ont_output:
    """
    This rule gzip the filtlong output reads into the correct location.
    """
    input:
        FASTQ = rules.trim_ont_workflow.output.FASTQ
    output:
        FASTQ = Path(config['working_dir']) / 'trimming' / 'ont' / '{}_SE.fastq.gz'.format(config['name'])
    params:
        dir_ = Path(config['working_dir']) / 'trimming' / 'ont'
    run:
        input_fastq = open(input.FASTQ, 'rb').read()
        with gzip.open(output.FASTQ, 'wb') as handle:
            handle.write(input_fastq)

rule unicycler:
    """
    Runs unicycler, which is a short-read first approach to assemble reads.
    """
    input:
        FQ_1P = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_1P.fastq.gz",
        FQ_2P = Path(config['working_dir']) / 'trimming' / 'illumina' / f"{config['name']}_2P.fastq.gz",
        FQ_SE = Path(config['working_dir']) / 'trimming' / 'ont' / '{}_SE.fastq.gz'.format(config['name'])
    output:
        FASTA = Path(config['working_dir']) / 'unicycler' / 'assembly.fasta',
        INFORMS = Path(config['working_dir']) / 'unicycler' / 'commands.io'
    params:
        working_dir = Path(config['working_dir'])
    threads: 12
    run:
        from camel.app.tools.unicycler.unicycler import Unicycler
        unicycler_assembly = Unicycler(camel)
        unicycler_assembly.add_input_files({'FASTQ_PE': [ToolIOFile(Path(input.FQ_1P)), ToolIOFile(Path(input.FQ_2P))],
                                            'FASTQ_SE': [ToolIOFile(Path(input.FQ_SE))]})
        unicycler_assembly.update_parameters(output_dir='unicycler', threads=threads)
        step = Step(str(rule), unicycler_assembly, camel, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_object(unicycler_assembly.informs, Path(output.INFORMS))

rule report_quast:
    """
    Rule to generate the QUAST report table.
    """
    input:
        report_quast = [Path(config['working_dir']) / 'qc' / f'{name}' / 'quast' / 'informs.io' for name in assembly_steps]
        # report_quast_combined = Path(config['working_dir']) / 'qc' / 'quast_combined' / 'informs.io'
    output:
        INFORMS = Path(config['working_dir'] / 'report' / 'quast.io')
    run:
        quast_table = []
        for i in range(len(input.report_quast)):
            quast_informs = SnakemakeUtils.load_object(Path(input.report_quast[i]))
            quast_table.extend([
                (assembly_steps[i], '{:,}'.format(int(quast_informs['contig']['N50'])),
                 '{:,}'.format(int(quast_informs['contig']['# contigs (>= 1000 bp)'])),
                 '{:,}'.format(int(quast_informs['genome']['Total length']))),
            ])
        SnakemakeUtils.dump_object(quast_table, Path(output.INFORMS))

rule report_short_variant_calling:
    """
    Rule to generate the report table for freebayes/clair3.
    """
    input:
        report_freebayes = [Path(config['working_dir']) / 'qc' / f'{name}' / 'freebayes' / 'informs.io' for name in assembly_steps],
        report_clair3 = [Path(config['working_dir']) / 'qc' / f'{name}' / 'clair3_output' / 'informs.io' for name in assembly_steps]
    output:
        INFORMS = Path(config['working_dir'] / 'report' / 'variant_calling.io')
    run:
        vc_table = []
        for i in range(len(input.report_freebayes)):
            freebayes_informs = SnakemakeUtils.load_object(Path(input.report_freebayes[i]))
            clair3_informs = SnakemakeUtils.load_object(Path(input.report_clair3[i]))
            vc_table.extend([
                (assembly_steps[i], '{:,}'.format(int(freebayes_informs['number_of_variants'])),
                 '{:,}'.format(int(len([t for t in freebayes_informs['type_of_variants'] if
                                        freebayes_informs['type_of_variants'] != [] and t.var_type == 'indel']))),
                 '{:,}'.format(int(len([t for t in freebayes_informs['type_of_variants'] if
                                        freebayes_informs['type_of_variants'] != [] and t.var_type == 'snp']))),
                 '{:,}'.format(int(clair3_informs['number_of_variants'])),
                 '{:,}'.format(int(len([t for t in clair3_informs['type_of_variants'] if
                                        clair3_informs['type_of_variants'] != [] and t.var_type == 'indel']))),
                 '{:,}'.format(int(len([t for t in clair3_informs['type_of_variants'] if
                                        clair3_informs['type_of_variants'] != [] and t.var_type == 'snp'])))
                 )
            ])
        SnakemakeUtils.dump_object(vc_table, Path(output.INFORMS))

rule report_long_variant_calling:
    """
    Rule to generate the report table for sniffles.
    """
    input:
        report_sniffles = [Path(config['working_dir']) / 'qc' / f'{name}' / 'sniffles' / 'informs.io' for name in assembly_steps]
    output:
        INFORMS = Path(config['working_dir'] / 'report' / 'sniffles.io')
    run:
        sniffles_table = []
        for i in range(len(input.report_sniffles)):
            sniffles_informs = SnakemakeUtils.load_object(Path(input.report_sniffles[i]))
            sniffles_table.extend([
                (assembly_steps[i], '{:,}'.format(int(sniffles_informs['number_of_variants'])),
                 '{:,}'.format(int(len([t for t in sniffles_informs['type_of_variants'] if
                                        sniffles_informs['type_of_variants'] != [] and t.var_type == 'indel']))),
                 '{:,}'.format(int(len([t for t in sniffles_informs['type_of_variants'] if
                                        sniffles_informs['type_of_variants'] != [] and t.var_type == 'snp'])))
                 )
            ])
        SnakemakeUtils.dump_object(sniffles_table, Path(output.INFORMS))

rule report_mappingstats:
    """
    Rule to generate the report table for the mapping statistics.
    """
    input:
        report_longreads = [Path(config['working_dir']) / 'qc' / f'{name}' / 'read_mapping' / 'ont' / 'flagstat-longreads.io'  for name in assembly_steps],
        report_shortreads = [Path(config['working_dir']) / 'qc' / f'{name}' / 'read_mapping' / 'illumina' / 'flagstat.io' for name in assembly_steps],
        report_depth_shortreads = [Path(config['working_dir']) / 'qc' / f'{name}' / 'read_mapping' / 'illumina' / 'samtools-depth.io' for name in assembly_steps],
        report_depth_longreads = [Path(config['working_dir']) / 'qc' / f'{name}' / 'read_mapping' / 'ont' / 'samtools-depth-long.io' for name in assembly_steps]
    output:
        INFORMS = Path(config['working_dir'] / 'report' / 'mapping.io')
    run:
        mapping_table = []
        for i in range(len(input.report_longreads)):
            SR_informs = SnakemakeUtils.load_object(Path(input.report_shortreads[i]))
            LR_informs = SnakemakeUtils.load_object(Path(input.report_longreads[i]))
            SR_depth = SnakemakeUtils.load_object(Path(input.report_depth_shortreads[i]))
            LR_depth = SnakemakeUtils.load_object(Path(input.report_depth_longreads[i]))
            mapping_rate_SR = SR_informs['mapped'][0] / SR_informs['total'][0] * 100
            mapping_rate_LR = LR_informs['mapped'][0] / LR_informs['total'][0] * 100
            median_depth_SR = int(SR_depth['median_depth'])
            median_depth_LR = int(LR_depth['median_depth'])
            print([mapping_rate_SR, mapping_rate_LR, median_depth_SR, median_depth_LR])
            mapping_table.extend([
                (assembly_steps[i], '{:.2f}%'.format(mapping_rate_SR), '{:,}'.format(median_depth_SR),
                 '{:.2f}%'.format(mapping_rate_LR), '{:,}'.format(median_depth_LR)
                 )
            ])
        SnakemakeUtils.dump_object(mapping_table,Path(output.INFORMS))

rule report_command_section:
    """
    Creates a report section with the commands used in the pipeline. 
    """
    input:
        unicycler_commands = Path(config['working_dir']) / 'unicycler' / 'commands.io',
        flye_commands = Path(config['working_dir']) / 'assembly_flye' / 'flye' / 'commands.io',
        medaka_consensus_commands = Path(config['working_dir']) / 'medaka' / 'commands-consensus.io',
        medaka_stitch_commands = Path(config['working_dir']) / 'medaka' / 'commands-stitch.io',
        polypolish_commands = Path(config['working_dir']) / 'polishing' / 'polypolish'  / 'polypolish.io',
        polca_commands = Path(config['working_dir']) / 'polishing' / 'polca' / 'polca.io',
        quast_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'quast' / 'commands.io' for name in assembly_steps],
        # quast_combined_commands = Path(config['working_dir']) / 'qc' / 'quast_combined' / 'commands.io',
        bwa_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'read_mapping' / 'illumina' / 'commands.io' for name in assembly_steps],
        freebayes_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'freebayes' / 'commands.io' for name in assembly_steps],
        sniffles_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'sniffles' / 'commands.io' for name in assembly_steps],
        clair3_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'clair3_output' / 'commands.io' for name in assembly_steps],
        ale_report_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'ale_illumina' / 'commands-report.io' for name in assembly_steps],
        ale_wiggle_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'ale_illumina' / 'commands-wiggle.io' for name in assembly_steps]
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-commands.io'
    params:
        working_dir = config['working_dir']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        informs = []
        for content in [SnakemakeUtils.load_object(Path(io)) for io in input]:
            if type(content) is dict:
                informs.append(content)
            elif type(content) is list:
                informs.extend(content)
        section = SnakePipelineUtils.create_commands_section(informs, params.working_dir)
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule report_create_sections:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        informs_quast = rules.report_quast.output.INFORMS,
        informs_vc = rules.report_short_variant_calling.output.INFORMS,
        informs_sniffles = rules.report_long_variant_calling.output.INFORMS,
        informs_mapping = rules.report_mappingstats.output.INFORMS,
        report_commands = rules.report_command_section.output.HTML,
        informs_trimming_illumina = rules.trim_illumina.output.HTML,
        informs_trimming_ont = rules.trim_ont_workflow.output.HTML
    output:
        HTML = Path(config['working_dir']) / Path(config['output_html'])
    params:
        sample_name = config['sample_name'],
        fastq_input = config['input']['illumina_name'],
        fastq_se_input = config['input']['ont_name'],
        pipeline = config['pipeline'],
        working_dir = Path(config['working_dir'])
    run:
        from camel.scripts.hybridassemblypipeline.reporter.hybridassemblyreporter import HybridAssemblyReporter
        hybridassembly_reporter = HybridAssemblyReporter(camel, params.working_dir)
        hybridassembly_reporter.add_input_informs({'quast': SnakemakeUtils.load_object(Path(input.informs_quast)),
                                                   'vc': SnakemakeUtils.load_object(Path(input.informs_vc)),
                                                   'sniffles': SnakemakeUtils.load_object(Path(input.informs_sniffles)),
                                                   'mapping': SnakemakeUtils.load_object(Path(input.informs_mapping)),
                                                   'commands': SnakemakeUtils.load_object(Path(input.report_commands)),
                                                   'trimming_illumina': SnakemakeUtils.load_object(Path(input.informs_trimming_illumina)),
                                                   'trimming_ont': SnakemakeUtils.load_object(Path(input.informs_trimming_ont)),
                                                   'sample_name': params.sample_name,
                                                   'fastq_input': params.fastq_input,
                                                   'fastq_se_input': params.fastq_se_input,
                                                   'pipeline': params.pipeline
                                                 })
        step = Step(str(rule), hybridassembly_reporter, camel, params.working_dir, config)
        step.run_step()
