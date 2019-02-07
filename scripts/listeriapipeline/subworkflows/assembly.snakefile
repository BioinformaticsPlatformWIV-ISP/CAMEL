ASSEMBLY_WORKING_DIR = os.path.join(__WORKING_DIR, 'assembly')
ASSEMBLY_REPORT = os.path.join(ASSEMBLY_WORKING_DIR, 'report-html.io')
ASSEMBLY_SUMMARY = os.path.join(ASSEMBLY_WORKING_DIR, 'report-summary.tsv')
FASTA_ASSEMBLY_RAW = os.path.join(ASSEMBLY_WORKING_DIR, 'fasta_raw.io')
FASTA_ASSEMBLY = os.path.join(ASSEMBLY_WORKING_DIR, 'fasta.io')

rule velvet_optimiser:
    """
    De-novo assembly using VelvetOptimiser.
    """
    input:
        FASTQ_PE = TRIMMED_READS_PE,
        FASTQ_SE_FORWARD = TRIMMED_READS_SE_FORWARD,
        FASTQ_SE_REVERSE = TRIMMED_READS_SE_REVERSE
    output:
        FASTA_Contig = os.path.join(ASSEMBLY_WORKING_DIR, 'velvet', 'fasta.io')
    params:
        running_dir = os.path.join(ASSEMBLY_WORKING_DIR, 'velvet'),
        fast = config.get('assembly_fast', False)
    threads:
        6
    run:
        from camel.app.tools.velvetoptimiser.velvetoptimiser import VelvetOptimiser
        velvet_optimiser = VelvetOptimiser(camel)
        SnakemakeUtils.add_pickle_inputs(velvet_optimiser, input)
        step = Step(rule, velvet_optimiser, camel, params.running_dir, config)
        velvet_optimiser.update_parameters(threads=threads)
        if params.fast:
            velvet_optimiser.update_parameters(hash_start=51, hash_end=51)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(velvet_optimiser, output)

rule spades:
    """
    De-novo assembly using SPAdes.
    """
    input:
        FASTQ_PE = TRIMMED_READS_PE,
        FASTQ_SE_FORWARD = TRIMMED_READS_SE_FORWARD,
        FASTQ_SE_REVERSE = TRIMMED_READS_SE_REVERSE
    output:
        FASTA_Contig = os.path.join(ASSEMBLY_WORKING_DIR, 'spades', 'fasta.io')
    params:
        running_dir = os.path.join(ASSEMBLY_WORKING_DIR, 'spades'),
        fast = config.get('assembly_fast', False)
    threads:
        6
    run:
        from camel.app.tools.spades.spades import SPAdes
        spades = SPAdes(camel)
        spades.add_input_files({
            'FASTQ_PE_1': SnakemakeUtils.load_object(input.FASTQ_PE),
            'FASTQ_PE-S_1': SnakemakeUtils.load_object(input.FASTQ_SE_FORWARD) +
            SnakemakeUtils.load_object(input.FASTQ_SE_REVERSE)
        })
        step = Step(rule, spades, camel, params.running_dir, config)
        spades.update_parameters(threads=threads, cov_cutoff='off')
        if params.fast:
            spades.update_parameters(only_assembly=True, kmers='55', careful=False)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spades, output)

rule select_assembly:
    """
    This rule selects the right assembler based on the config file.
    """
    input:
        os.path.join(ASSEMBLY_WORKING_DIR, 'spades', 'fasta.io') if config.get('assembler') == 'SPAdes' else [],
        os.path.join(ASSEMBLY_WORKING_DIR, 'velvet', 'fasta.io') if config.get('assembler') == 'VelvetOptimiser' else []
    output:
        FASTA_ASSEMBLY_RAW
    run:
        import shutil
        shutil.copyfile(input[0], output[0])

rule contig_filtering:
    input:
        FASTA_contig = FASTA_ASSEMBLY_RAW
    output:
        FASTA_ASSEMBLY
    params:
        running_dir = ASSEMBLY_WORKING_DIR
    run:
        from camel.app.tools.pipelines.assembly.contigqcfilter import ContigQCFilter
        ctgfilter = ContigQCFilter(camel)
        SnakemakeUtils.add_pickle_inputs(ctgfilter, input)
        step = Step(rule, ctgfilter, camel, params.running_dir, config)
        # ctgfilter.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(ctgfilter, output)
        
rule quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = FASTA_ASSEMBLY
    output:
        TSV = os.path.join(ASSEMBLY_WORKING_DIR, 'quast', 'tsv.io')
    params:
        running_dir = os.path.join(ASSEMBLY_WORKING_DIR, 'quast')
    threads:
        8
    run:
        from camel.app.tools.quast.quast import Quast
        quast = Quast(camel)
        SnakemakeUtils.add_pickle_inputs(quast, input)
        step = Step(rule, quast, camel, params.running_dir, config)
        quast.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast, output)

rule quast_inform_extractor:
    """
    Extracts the information from the QUAST output file.
    """
    input:
        TSV = os.path.join(ASSEMBLY_WORKING_DIR, 'quast', 'tsv.io')
    output:
        INFORMS = os.path.join(ASSEMBLY_WORKING_DIR, 'quast', 'informs.io')
    params:
        running_dir = os.path.join(ASSEMBLY_WORKING_DIR, 'quast')
    run:
        from camel.app.tools.quast.quastinformextractor import QuastInformExtractor
        quast_inform_extractor = QuastInformExtractor(camel)
        SnakemakeUtils.add_pickle_inputs(quast_inform_extractor, input)
        step = Step(rule, quast_inform_extractor, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast_inform_extractor, output)

rule report_assembly:
    """
    Creates the HTML report for the assembly.
    """
    input:
        FASTA_Contig = FASTA_ASSEMBLY,
        INFORMS_quast = os.path.join(ASSEMBLY_WORKING_DIR, 'quast', 'informs.io')
    output:
        VAL_HTML = ASSEMBLY_REPORT
    params:
        running_dir = os.path.join(ASSEMBLY_WORKING_DIR),
        sample_name = config['sample_name'],
        assembler = config['assembler'],
        output_dir = config['output_dir']
    run:
        from camel.app.tools.pipelines.assembly.htmlreporterassembly import HtmlReporterAssembly
        reporter = HtmlReporterAssembly(camel)
        reporter.add_input_files(
            {'SAMPLE_NAME': [ToolIOValue(params.sample_name)],
             'ASSEMBLER': [ToolIOValue(params.assembler)]})
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        reporter.tool_outputs['VAL_HTML'][0].value.copy_files(params.output_dir)
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule summary_assembly:
    """
    Creates a tabular summary for the assembly.
    """
    input:
        INFORMS_quast = os.path.join(ASSEMBLY_WORKING_DIR, 'quast', 'informs.io')
    output:
        ASSEMBLY_SUMMARY
    params:
        running_dir = os.path.join(ASSEMBLY_WORKING_DIR),
        sample_name = config['sample_name'],
        assembler = config['assembler'],
        output_dir = config['output_dir']
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
