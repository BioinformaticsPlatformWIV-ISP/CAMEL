from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import quast, assembly


rule quast_quast:
    """
    Runs quast on the assembly.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.get_fasta(config),
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        TSV = Path(config['working_dir']) / 'quast' / 'output' / 'tsv.io',
        HTML = Path(config['working_dir']) / 'quast' / 'output' / 'html.io',
        DIR = Path(config['working_dir']) / 'quast' / 'output' / 'dir.io',
        INFORMS = Path(config['working_dir']) / 'quast' / 'output' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'quast' / 'output',
        input_type = config['input_type'],
        fasta = config.get('quast', {}).get('ref', {}).get('fasta'),
        gff = config.get('quast', {}).get('ref', {}).get('gff3')
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.tools.quast.quast import Quast

        # Create output directory
        dir_out = Path(params.running_dir) / 'quast_out'
        dir_out.mkdir(exist_ok=True, parents=True)

        # Create tool
        quast_ = Quast(Camel.get_instance())

        # Add input
        SnakemakeUtils.add_pickle_inputs(quast_, input, excluded_keys=['IO'])
        fq_dict = SnakemakeUtils.load_object(Path(input.IO))
        if params.input_type in ('illumina', 'hybrid'):
            quast_.add_input_files({'FASTQ_PE': fq_dict['PE']})
        if params.input_type in ('ont', 'hybrid'):
            quast_.add_input_files({'FASTQ_nanopore': fq_dict['SE']})

        # Add reference genome files (if available)
        if (params.fasta is not None) and (params.gff is not None):
            quast_.add_input_files({
                'FASTA_Ref': [ToolIOFile(Path(params.fasta))],
                'GFF3_Ref': [ToolIOFile(Path(params.gff))],
            })
        else:
            logger.warning(f'No reference genome provided, skipping analysis for QUAST')

        # Run tool
        quast_.update_parameters(conserved_genes_finding=False)
        step = Step(str(rule), quast_, Camel.get_instance(), dir_out)
        step.run_step()

        # Collect output
        SnakemakeUtils.dump_tool_outputs(quast_, output, ignore_missing_output=True)
        SnakemakeUtils.dump_object([ToolIODirectory(dir_out)], Path(output.DIR))

rule quast_busco:
    """
    Runs BUSCO on the assembly to check completeness.
    BUSCO is ran outside of QUAST because of dependency issues with the BUSCO installation bundles with QUAST.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.get_fasta(config)
    output:
        TXT = Path(config['working_dir']) / 'quast' / 'busco' / 'txt.io',
        INFORMS = Path(config['working_dir']) / 'quast' / 'busco' / 'informs.io'
    params:
        dir_ = Path(config['working_dir']) / 'quast' / 'busco',
        lineage_dataset = 'bacteria_odb10'
    threads: 8
    run:
        from camel.app.tools.busco.busco import Busco
        busco = Busco(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(busco, input)
        step = Step(str(rule), busco, Camel.get_instance(), params.dir_)
        busco.update_parameters(lineage_dataset=params.lineage_dataset, threads=str(threads))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(busco, output)

rule quast_report:
    """
    Creates a report for the quast workflow.
    """
    input:
        TSV = rules.quast_quast.output.TSV,
        HTML = rules.quast_quast.output.HTML,
        FASTA =  Path(config['working_dir']) / assembly.get_fasta(config),
        DIR = rules.quast_quast.output.DIR,
        INFORMS_quast = rules.quast_quast.output.INFORMS,
        INFORMS_assembler = Path(config['working_dir']) / assembly.get_command_informs(config),
        INFORMS_busco = rules.quast_busco.output.INFORMS
    output:
        HTML = Path(config['working_dir']) / 'quast' / 'report' / 'html.io'
    params:
        running_dir = Path(config['working_dir']) / 'quast' / 'report',
        name = config['sample_name']
    run:
        from camel.app.tools.quast.quastreporter import QuastReporter
        reporter = QuastReporter(Camel.get_instance())
        reporter.update_parameters(name=params.name)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, Camel.get_instance(), params.running_dir)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule quast_create_summary_out:
    """
    Creates the tabular summary output for QUAST.
    """
    input:
        INFORMS = Path(config['working_dir'], assembly.get_command_informs(config)),
        TSV = rules.quast_quast.output.TSV
    output:
        TSV = Path(config['working_dir']) / quast.OUTPUT_QUAST_SUMMARY
    run:
        keys_kept = [
            {'key': '# contigs', 'name': 'nb_contigs'},
            {'key': 'Total length', 'name': 'total_length'},
            {'key': 'Reference length', 'name': 'total_length_ref'},
            {'key': 'N50', 'name': 'n50'},
            {'key': 'Genome fraction (%)', 'name': 'genome_fraction'},
            {'key': 'Duplication ratio', 'name': 'dupl_ratio'},
            {'key': 'tool_version', 'name': 'tool_version'},
        ]

        # Parse QUAST report
        path_tsv = SnakemakeUtils.load_object(Path(input.TSV))[0].path
        data_quast = {}
        with path_tsv.open() as handle:
            for line in handle.readlines():
                key, value = line.strip().split('\t')
                data_quast[key] = value

        # Add assembler version
        spades_informs = SnakemakeUtils.load_object(Path(input.INFORMS))
        data_quast['tool_version'] = spades_informs['_name']

        # Create TSV output
        with open(output.TSV, 'w') as handle:
            for row in keys_kept:
                handle.write(f"assembly_{row['name']}\t{data_quast.get(row['key'], '-')}")
                handle.write('\n')
