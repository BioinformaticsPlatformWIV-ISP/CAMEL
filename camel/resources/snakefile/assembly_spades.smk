import os

from camel.app.camel import Camel
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile.assembly_spades import OUTPUT_ASSEMBLY_SUMMARY
from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_READS_PE, OUTPUT_READ_TRIMMING_READS_SE_FWD, \
    OUTPUT_READ_TRIMMING_READS_SE_REV

camel = Camel.get_instance()


rule Assembly_spades:
    """
    De-novo assembly using SPAdes.
    """
    input:
        FASTQ_PE=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_PE),
        FASTQ_SE_FORWARD=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_SE_FWD),
        FASTQ_SE_REVERSE=os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_READS_SE_REV)
    output:
        FASTA_Contig=os.path.join(config['working_dir'], 'assembly_spades', 'spades', 'fasta.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'assembly_spades', 'spades'),
        kmers=config['assembly'].get('kmers') if 'assembly' in config else None
    threads: 8
    run:
        from camel.app.tools.spades.spades import SPAdes
        spades = SPAdes(camel)
        spades.add_input_files({
            'FASTQ_PE_1': SnakemakeUtils.load_object(input.FASTQ_PE),
            'FASTQ_PE-S_1': SnakemakeUtils.load_object(input.FASTQ_SE_FORWARD) +
                            SnakemakeUtils.load_object(input.FASTQ_SE_REVERSE)
        })
        step = Step(rule, spades, camel, params.running_dir, config)
        if params.kmers is not None:
            spades.update_parameters(kmers=params.kmers)
        spades.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spades, output)

rule Assembly_quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA=os.path.join(config['working_dir'], 'assembly_spades', 'spades', 'fasta.io')
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
        FASTA_Contig=os.path.join(config['working_dir'], 'assembly_spades', 'spades', 'fasta.io'),
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
    Dumps the summary information from the read trimming pipeline.
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
