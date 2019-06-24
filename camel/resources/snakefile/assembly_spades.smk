import os

from camel.app.camel import Camel
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.assembly_spades import OUTPUT_ASSEMBLY_SUMMARY, OUTPUT_ASSEMBLY_INFORMS
from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_READS_PE, OUTPUT_READ_TRIMMING_READS_SE_FWD, \
    OUTPUT_READ_TRIMMING_READS_SE_REV
from camel.resources.snakefile.read_trimming_iontorrent import OUTPUT_TRIMMING_IT_READS

camel = Camel.get_instance()


rule Assembly_select_fastq_input:
    """
    This rule is used to select the assembly input based on the read type input specified in the config.
    'illumina' (default): Trimmed PE reads + orphaned SE reads (optional)
    'iontorrent': Trimmed SE reads
    """
    input:
        ILLUMINA_FASTQ_PE=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_PE) if config.get('read_type', 'illumina') == 'illumina' else [],
        ILLUMINA_FASTQ_SE_FWD=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_SE_FWD) if config.get('read_type', 'illumina') == 'illumina' else [],
        ILLUMINA_FASTQ_SE_REV=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_SE_REV) if config.get('read_type', 'illumina') == 'illumina' else [],
        IONTORRENT_FASTQ_SE=os.path.join(config['working_dir'], OUTPUT_TRIMMING_IT_READS) if config.get('read_type', 'illumina') == 'iontorrent' else []
    output:
        os.path.join(config['working_dir'], 'assembly_spades', 'spades', 'input.io')
    params:
        read_type=config.get('read_type', 'illumina')
    run:
        output_dict = {}
        if params.read_type == 'illumina':
            output_dict = {'FASTQ_PE_1': SnakemakeUtils.load_object(input.ILLUMINA_FASTQ_PE)}
            se_reads = SnakemakeUtils.load_object(input.ILLUMINA_FASTQ_SE_FWD) + \
                       SnakemakeUtils.load_object(input.ILLUMINA_FASTQ_SE_REV)
            if len(se_reads) > 0:
                output_dict['FASTQ_SE_1'] = se_reads
        else:
            tmp = SnakemakeUtils.load_object(input.IONTORRENT_FASTQ_SE)
            output_dict = {'FASTQ_SE_1': SnakemakeUtils.load_object(input.IONTORRENT_FASTQ_SE)}
        SnakemakeUtils.dump_object(output_dict, output[0])

rule Assembly_spades:
    """
    De-novo assembly using SPAdes.
    """
    input:
        INPUT_DICT=os.path.join(config['working_dir'], 'assembly_spades', 'spades', 'input.io')
    output:
        FASTA_Contig=os.path.join(config['working_dir'], 'assembly_spades', 'spades', 'fasta.io'),
        INFORMS=os.path.join(config['working_dir'], OUTPUT_ASSEMBLY_INFORMS)
    params:
        running_dir=os.path.join(config['working_dir'], 'assembly_spades', 'spades'),
        spades_options=config.get('assembly', {}).get('spades', {})
    threads: 8
    priority: 1
    run:
        from camel.app.tools.spades.spades import SPAdes
        spades = SPAdes(camel)
        spades.add_input_files(SnakemakeUtils.load_object(input.INPUT_DICT))
        step = Step(rule, spades, camel, params.running_dir, config)
        spades.update_parameters(**params.spades_options)
        spades.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spades, output)

rule Assembly_filter_small_contigs:
    """
    Filters out the small contigs.
    """
    input:
        FASTA = os.path.join(config['working_dir'], 'assembly_spades', 'spades', 'fasta.io')
    output:
        FASTA = os.path.join(config['working_dir'], 'assembly_spades', 'filtering', 'fasta.io'),
        INFORMS = os.path.join(config['working_dir'], 'assembly_spades', 'filtering', 'informs.io')
    params:
        running_dir = os.path.join(config['working_dir'], 'assembly_spades', 'filtering'),
        min_contig_length = config['assembly'].get('min_contig_length', 0) if 'assembly' in config else 0
    run:
        from camel.app.tools.seqtk.seqtkseq import SeqtkSeq
        seqtk = SeqtkSeq(camel)
        SnakemakeUtils.add_pickle_inputs(seqtk, input)
        step = Step(rule, seqtk, camel, params.running_dir, config)
        seqtk.update_parameters(output_filename='assembly_filtered.fasta', min_length=params.min_contig_length)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqtk, output)

rule Assembly_quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = os.path.join(config['working_dir'], 'assembly_spades', 'filtering', 'fasta.io')
    output:
        TSV=os.path.join(config['working_dir'], 'assembly_spades', 'quast', 'tsv.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'assembly_spades', 'quast')
    run:
        from camel.app.tools.quast.quast import Quast
        quast = Quast(camel)
        SnakemakeUtils.add_pickle_inputs(quast, input)
        step = Step(rule, quast, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast, output)

rule Assembly_quast_inform_extractor:
    """
    Extracts the information from the QUAST output file.
    """
    input:
        TSV=os.path.join(config['working_dir'], 'assembly_spades', 'quast', 'tsv.io')
    output:
        INFORMS=os.path.join(config['working_dir'], 'assembly_spades', 'quast', 'informs.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'assembly_spades', 'quast')
    run:
        from camel.app.tools.quast.quastinformextractor import QuastInformExtractor
        quast_inform_extractor = QuastInformExtractor(camel)
        SnakemakeUtils.add_pickle_inputs(quast_inform_extractor, input)
        step = Step(rule, quast_inform_extractor, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast_inform_extractor, output)

rule Assembly_report:
    """
    Creates the HTML report for the assembly.
    """
    input:
        FASTA_Contig = os.path.join(config['working_dir'], 'assembly_spades', 'filtering', 'fasta.io'),
        INFORMS_spades=os.path.join(config['working_dir'], 'assembly_spades', 'spades', 'informs.io'),
        INFORMS_quast=os.path.join(config['working_dir'], 'assembly_spades', 'quast', 'informs.io')
    output:
        VAL_HTML=os.path.join(config['working_dir'], 'assembly_spades', 'report', 'html.io')
    params:
        running_dir = os.path.join(config['working_dir'], 'assembly_spades', 'report'),
        sample_name=config['sample_name']
    run:
        from camel.app.tools.pipelines.assembly.htmlreporterassembly import HtmlReporterAssembly
        reporter = HtmlReporterAssembly(camel)
        reporter.add_input_files({'SAMPLE_NAME': [ToolIOValue(params.sample_name)],
                                  'ASSEMBLER': [ToolIOValue('SPAdes')]})
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule Assembly_dump_summary_info:
    """
    Dumps the summary information from the assembly pipeline.
    """
    input:
        INFORMS_quast=os.path.join(config['working_dir'], 'assembly_spades', 'quast', 'informs.io')
    output:
        os.path.join(config['working_dir'], OUTPUT_ASSEMBLY_SUMMARY)
    params:
        running_dir=os.path.join(config['working_dir'], 'assembly_spades', 'summary')
    run:
        quast_informs = SnakemakeUtils.load_object(input.INFORMS_quast)
        summary_data = [
            ('n50', quast_informs['contig']['N50']),
            ('nb_contigs', quast_informs['contig']['# contigs']),
            ('nb_contigs_lt_1000', quast_informs['contig']['# contigs (>= 1000 bp)']),
            ('total_length', quast_informs['genome']['Total length'])
        ]
        with open(output[0], 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')
