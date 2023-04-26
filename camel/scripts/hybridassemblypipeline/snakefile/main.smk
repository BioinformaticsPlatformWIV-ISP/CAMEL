import gzip
import pickle
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

#########
# Rules #
#########

rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        Path(config['working_dir']) / config['output_html'],
        Path(config['working_dir'] / config['output']),
        Path(config['working_dir']) / 'report' / 'commands.txt'

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
        HTML = Path(config['working_dir']) / 'trimming' / 'illumina' / 'fastqc_pre.html'
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
        wrapper.output.fastq_reports_pre[0].path.rename(Path(output.HTML))

rule trim_ont:
    """
    This rule trims the ONT reads using filtlong.
    """
    input:
        FASTQ = config['input']['ont']
    output:
        FASTQ = Path(config['working_dir']) / 'trimming' / 'ont' / 'fastq.io'
    params:
        dir_ = Path(config['working_dir']) / 'trimming' / 'ont',
        filtlong_options= config.get('filtlong',{})
    threads: 4
    run:
        from camel.app.tools.filtlong.filtlong import Filtlong
        filtlong = Filtlong(camel)
        filtlong.add_input_files({'FASTQ': [ToolIOFile(Path(input.FASTQ))]})
        filtlong.update_parameters(**params.filtlong_options)
        step = Step(str(rule), filtlong, camel, params.dir_)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(filtlong,output)

rule set_trimming_ont_output:
    """
    This rule gzip the filtlong output reads into the correct location.
    """
    input:
        FASTQ = rules.trim_ont.output.FASTQ
    output:
        FASTQ = Path(config['working_dir']) / 'trimming' / 'ont' / '{}_SE.fastq.gz'.format(config['name'])
    params:
        dir_ = Path(config['working_dir']) / 'trimming' / 'ont'
    run:
        input_fastq = open(SnakemakeUtils.load_object(Path(input.FASTQ))[0].path, 'rb').read()
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
        dir_ = Path(config['working_dir'])
    threads: 4
    run:
        from camel.app.tools.unicycler.unicycler import Unicycler
        unicycler_assembly = Unicycler(camel)
        unicycler_assembly.add_input_files({'FASTQ_PE':[ToolIOFile(Path(input.FQ_1P)), ToolIOFile(Path(input.FQ_2P))],
                                            'FASTQ_SE':[ToolIOFile(Path(input.FQ_SE))]})
        unicycler_assembly.update_parameters(output_dir = 'unicycler', threads=threads)
        step = Step(str(rule), unicycler_assembly, camel, params.dir_, config)
        step.run_step()
        with open(output.INFORMS, 'wb') as handle:
            pickle.dump(unicycler_assembly.informs, handle)

rule report_quast:
    """
    Rule to generate the QUAST report table.
    """
    input:
        report_quast = [Path(config['working_dir']) / 'qc' / f'{name}' / 'quast' / 'informs.io' for name in ['medaka', 'polca', 'polypolish', 'unicycler']]
    output:
        INFORMS = Path(config['working_dir'] / 'report' / 'quast.io')
    run:
        names = ['medaka', 'polca', 'polypolish', 'unicycler']
        quast_table = []
        for i in range(len(input.report_quast)):
            quast_informs = SnakemakeUtils.load_object(Path(input.report_quast[i]))
            quast_table.extend([
                (names[i], '{:,}'.format(int(quast_informs['contig']['N50'])),
                 '{:,}'.format(int(quast_informs['contig']['# contigs (>= 1000 bp)'])),
                 '{:,}'.format(int(quast_informs['genome']['Total length']))),
            ])
        with open(output.INFORMS,'wb') as handle:
            pickle.dump(quast_table,handle)

rule report_short_variant_calling:
    """
    Rule to generate the report table for freebayes/clair3.
    """
    input:
        report_freebayes = [Path(config['working_dir']) / 'qc' / f'{name}' / 'freebayes' / 'informs.io' for name in ['medaka', 'polca', 'polypolish', 'unicycler']],
        report_clair3 = [Path(config['working_dir']) / 'qc' / f'{name}' / 'clair3_output' / 'informs.io' for name in ['medaka', 'polca', 'polypolish', 'unicycler']]
    output:
        INFORMS = Path(config['working_dir'] / 'report' / 'variant_calling.io')
    run:
        names = ['medaka', 'polca', 'polypolish', 'unicycler']
        vc_table = []
        for i in range(len(names)):
            freebayes_informs = SnakemakeUtils.load_object(Path(input.report_freebayes[i]))
            clair3_informs = SnakemakeUtils.load_object(Path(input.report_clair3[i]))
            vc_table.extend([
                (names[i], '{:,}'.format(int(freebayes_informs['number_of_variants'])),
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
        with open(output.INFORMS,'wb') as handle:
            pickle.dump(vc_table,handle)

rule report_long_variant_calling:
    """
    Rule to generate the report table for sniffles.
    """
    input:
        report_sniffles = [Path(config['working_dir']) / 'qc' / f'{name}' / 'sniffles' / 'informs.io' for name in ['medaka', 'polca', 'polypolish', 'unicycler']]
    output:
        INFORMS = Path(config['working_dir'] / 'report' / 'sniffles.io')
    run:
        names = ['medaka', 'polca', 'polypolish', 'unicycler']
        sniffles_table = []
        for i in range(len(input.report_sniffles)):
            sniffles_informs = SnakemakeUtils.load_object(Path(input.report_sniffles[i]))
            sniffles_table.extend([
                (names[i], '{:,}'.format(int(sniffles_informs['number_of_variants'])),
                 '{:,}'.format(int(len([t for t in sniffles_informs['type_of_variants'] if
                                        sniffles_informs['type_of_variants'] != [] and t.var_type == 'indel']))),
                 '{:,}'.format(int(len([t for t in sniffles_informs['type_of_variants'] if
                                        sniffles_informs['type_of_variants'] != [] and t.var_type == 'snp'])))
                 )
            ])
        with open(output.INFORMS,'wb') as handle:
            pickle.dump(sniffles_table,handle)

rule report_mappingstats:
    input:
        report_longreads = [Path(config['working_dir']) / 'qc' / f'{name}' / 'read_mapping' / 'flagstat-longreads.io'  for name in ['medaka', 'polca', 'polypolish', 'unicycler']],
        report_shortreads = [Path(config['working_dir']) / 'qc' / f'{name}' / 'read_mapping' / 'flagstat.io' for name in ['medaka', 'polca', 'polypolish', 'unicycler']],
        report_depth_shortreads = [Path(config['working_dir']) / 'qc' / f'{name}' / 'read_mapping' / 'samtools-depth.io' for name in ['medaka', 'polca', 'polypolish', 'unicycler']],
        report_depth_longreads = [Path(config['working_dir']) / 'qc' / f'{name}' / 'read_mapping' / 'samtools-depth-long.io' for name in ['medaka', 'polca', 'polypolish', 'unicycler']]
    output:
        INFORMS = Path(config['working_dir'] / 'report' / 'mapping.io')
    run:
        names = ['medaka', 'polca', 'polypolish', 'unicycler']
        mapping_table = []
        for i in range(len(input.report_longreads)):
            SR_informs = SnakemakeUtils.load_object(Path(input.report_shortreads[i]))
            LR_informs = SnakemakeUtils.load_object(Path(input.report_longreads[i]))
            SR_depth = SnakemakeUtils.load_object(Path(input.report_depth_shortreads[i]))
            LR_depth = SnakemakeUtils.load_object(Path(input.report_depth_longreads[i]))
            mapping_rate_SR = SR_informs['mapped'][0] / SR_informs['total'][0]
            mapping_rate_LR = LR_informs['mapped'][0] / LR_informs['total'][0]
            median_depth_SR = SR_depth['median_depth']
            median_depth_LR = LR_depth['median_depth']
            mapping_table.extend([
                (names[i], '{:.2f}'.format(mapping_rate_SR), '{:,}'.format(median_depth_SR),
                 '{:.2f}'.format(mapping_rate_LR), '{:,}'.format(median_depth_LR)
                 )
            ])
            with open(output.INFORMS, 'wb') as handle:
                pickle.dump(mapping_table, handle)

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
        quast_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'quast' / 'commands.io' for name in
                        ['medaka', 'polca', 'polypolish', 'unicycler']],
        bwa_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'read_mapping' / 'commands.io' for name in
                      ['medaka', 'polca', 'polypolish', 'unicycler']],
        freebayes_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'freebayes' / 'commands.io' for name in
                      ['medaka', 'polca', 'polypolish', 'unicycler']],
        sniffles_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'sniffles' / 'commands.io' for name in
                      ['medaka', 'polca', 'polypolish', 'unicycler']],
        clair3_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'clair3_output' / 'commands.io' for name in
                      ['medaka', 'polca', 'polypolish', 'unicycler']],
        ale_report_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'ale_illumina' / 'commands-report.io' for name in
                      ['medaka', 'polca', 'polypolish', 'unicycler']],
        ale_wiggle_commands = [Path(config['working_dir']) / 'qc' / f'{name}' / 'ale_illumina' / 'commands-wiggle.io' for name in
                      ['medaka', 'polca', 'polypolish', 'unicycler']]
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

rule report_combine_all:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        informs_quast = rules.report_quast.output.INFORMS,
        informs_vc = rules.report_short_variant_calling.output.INFORMS,
        informs_sniffles = rules.report_long_variant_calling.output.INFORMS,
        informs_mapping = rules.report_mappingstats.output.INFORMS,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = Path(config['working_dir']) / config['output_html']
    params:
        sample_name = config['sample_name'],
        fastq_input = config['input']['illumina'],
        fastq_se_input = config['input']['ont'],
        output_dir = config['working_dir'],
        pipeline = config['pipeline']
    run:
        import datetime
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        # Add header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            params.sample_name, datetime.datetime.now(), '0.1',
            ', '.join([str(params.fastq_se_input)]+[str(entry) for entry in params.fastq_input]), []))

        report.add_header('QUAST statistics', 2)
        quast_table = pickle.load(open(input.informs_quast, 'rb'))
        report.add_table(quast_table, column_names=['Assembly step', 'N50', 'No of contigs', 'Total length'], table_attributes=[('class', 'information')])

        report.add_header('Variant calling statistics', 2)
        vc_table = pickle.load(open(input.informs_vc, 'rb'))
        report.add_table(vc_table, column_names=['Assembly step', 'Number of variants', 'Number of indels', 'Number of SNPs', 'Clair3 total variant', 'Clair3 indels', 'Clair3 SNPs'], table_attributes=[('class', 'information')])

        report.add_header('Sniffles statistics', 2)
        sniffles_table = pickle.load(open(input.informs_sniffles, 'rb'))
        report.add_table(sniffles_table, column_names=['Assembly step', 'Number of variants', 'Number of indels', 'Number of SNPs'], table_attributes=[('class', 'information')])

        report.add_header('Mapping statistics', 2)
        mapping_table = pickle.load(open(input.informs_mapping, 'rb'))
        report.add_table(mapping_table, column_names=['Assembly step', 'mapping rate (short reads)', 'median depth (short reads)', 'mapping rate (long reads)', 'median depth (long reads)'], table_attributes=[('class', 'information')])

        report.add_header('Commands', 2)
        commands_content = [('Commands', 'commands', [Path(input.report_commands)])]
        SnakePipelineUtils.add_report_content(report, commands_content)

        report.save()
