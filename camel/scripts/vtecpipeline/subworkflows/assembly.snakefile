rule velvet_optimiser:
    """
    De-novo assembly using VelvetOptimiser.
    """
    input:
        FASTQ_PE=os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-pe.io'),
        FASTQ_SE_FORWARD=os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-se-forward.io'),
        FASTQ_SE_REVERSE=os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-se-reverse.io'),
    output:
        FASTA_Contig=os.path.join(__WORKING_DIR, 'assembly_velvet', 'fasta.io')
    params:
        running_dir=os.path.join(__WORKING_DIR, 'assembly_velvet'),
        fast = config.get('assembly_fast', False)
    threads: 8
    run:
        from camel.app.tools.velvetoptimiser.velvetoptimiser import VelvetOptimiser
        velvet_optimiser = VelvetOptimiser(camel)
        SnakemakeUtils.add_pickle_inputs(velvet_optimiser, input)
        step = SnakeStep(rule, velvet_optimiser, camel, params.running_dir, config)
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
        FASTQ_PE=os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-pe.io'),
        FASTQ_SE_FORWARD=os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-se-forward.io'),
        FASTQ_SE_REVERSE=os.path.join(__WORKING_DIR, 'read_trimming', 'fastq-se-reverse.io'),
    output:
        FASTA_Contig=os.path.join(__WORKING_DIR, 'assembly_spades', 'fasta.io')
    params:
        running_dir=os.path.join(__WORKING_DIR, 'assembly_spades'),
        fast=config.get('assembly_fast', False)
    threads: 8
    run:
        from camel.app.tools.spades.spades import SPAdes
        spades = SPAdes(camel)
        spades.add_input_files({
            'FASTQ_PE_1': SnakemakeUtils.load_object(input.FASTQ_PE),
            'FASTQ_PE-S_1': SnakemakeUtils.load_object(input.FASTQ_SE_FORWARD) +
                            SnakemakeUtils.load_object(input.FASTQ_SE_REVERSE)
        })
        step = SnakeStep(rule, spades, camel, params.running_dir, config)
        spades.update_parameters(threads=threads)
        if params.fast:
            spades.update_parameters(only_assembly=True, kmers='55', careful=False)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(spades, output)

rule select_assembly:
    """
    This rule selects the right assembler based on the config file.
    """
    input:
        os.path.join(__WORKING_DIR, 'assembly_spades', 'fasta.io') if config.get('assembler') == 'SPAdes' else [],
        os.path.join(__WORKING_DIR, 'assembly_velvet', 'fasta.io') if config.get('assembler') == 'VelvetOptimiser' else []
    output:
        os.path.join(__WORKING_DIR, 'assembly', 'fasta.io')
    run:
        import shutil
        shutil.copyfile(input[0], output[0])

rule quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA=os.path.join(__WORKING_DIR, 'assembly', 'fasta.io')
    output:
        TSV=os.path.join(__WORKING_DIR, 'quast', 'tsv.io')
    params:
        running_dir=os.path.join(__WORKING_DIR, 'quast')
    run:
        from camel.app.tools.quast.quast import Quast
        quast = Quast(camel)
        SnakemakeUtils.add_pickle_inputs(quast, input)
        step = SnakeStep(rule, quast, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast, output)

rule quast_inform_extractor:
    """
    Extracts the information from the QUAST output file.
    """
    input:
        TSV=os.path.join(__WORKING_DIR, 'quast', 'tsv.io')
    output:
        INFORMS=os.path.join(__WORKING_DIR, 'quast', 'informs.io')
    params:
        running_dir=os.path.join(__WORKING_DIR, 'quast')
    run:
        from camel.app.tools.quast.quastinformextractor import QuastInformExtractor
        quast_inform_extractor = QuastInformExtractor(camel)
        SnakemakeUtils.add_pickle_inputs(quast_inform_extractor, input)
        step = SnakeStep(rule, quast_inform_extractor, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast_inform_extractor, output)

rule report_assembly:
    """
    Creates the HTML report for the assembly.
    """
    input:
        FASTA_Contig=os.path.join(__WORKING_DIR, 'assembly', 'fasta.io'),
        INFORMS_quast=os.path.join(__WORKING_DIR, 'quast', 'informs.io')
    output:
        VAL_HTML=os.path.join(__WORKING_DIR, 'report_assembly', 'html.io')
    params:
        running_dir = os.path.join(__WORKING_DIR, 'report_assembly'),
        output_dir=config['output_dir'],
        sample_name=config['sample_name'],
        assembler=config['assembler']
    run:
        from camel.app.tools.pipelines.assembly.htmlreporterassembly import HtmlReporterAssembly
        reporter = HtmlReporterAssembly(camel)
        reporter.add_input_files(
            {'SAMPLE_NAME': [ToolIOValue(params.sample_name)],
             'ASSEMBLER': [ToolIOValue(params.assembler)]})
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = SnakeStep(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        reporter.tool_outputs['VAL_HTML'][0].value.copy_files(params.output_dir)
        SnakemakeUtils.dump_tool_outputs(reporter, output)
