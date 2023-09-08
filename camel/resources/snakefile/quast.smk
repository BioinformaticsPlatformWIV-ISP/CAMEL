from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import assembly_spades, quast


rule quast_pickle_genome:
    """
    Creates an IO pickle for the reference genome.
    """
    input:
        FASTA = config.get('quast', {}).get('ref', {}).get('fasta'),
        GFF3 = config.get('quast', {}).get('ref', {}).get('gff3')
    output:
        FASTA = Path(config['working_dir']) / 'quast' / 'ref_genome' / 'fasta.io',
        GFF3 = Path(config['working_dir']) / 'quast' / 'ref_genome' / 'gff3.io'
    run:
        from camel.app.io.tooliofile import ToolIOFile
        SnakemakeUtils.dump_object([ToolIOFile(Path(input.FASTA))], Path(output.FASTA))
        SnakemakeUtils.dump_object([ToolIOFile(Path(input.GFF3))], Path(output.GFF3))

rule quast_quast:
    """
    Runs quast on the assembly.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA,
        FASTA_Ref = Path(config['working_dir']) / rules.quast_pickle_genome.output.FASTA,
        GFF3_Ref = Path(config['working_dir']) / rules.quast_pickle_genome.output.GFF3,
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        TSV = Path(config['working_dir']) / 'quast' / 'output' / 'tsv.io',
        HTML = Path(config['working_dir']) / 'quast' / 'output' / 'html.io',
        DIR = Path(config['working_dir']) / 'quast' / 'output' / 'dir.io',
        INFORMS = Path(config['working_dir']) / 'quast' / 'output' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'quast' / 'output',
        read_type = config.get('read_type','illumina')
    run:
        from camel.app.tools.quast.quast import Quast

        # Create output directory
        dir_out = Path(params.running_dir) / 'quast_out'
        dir_out.mkdir(exist_ok=True, parents=True)

        # Run tool
        quast_ = Quast(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(quast_, input, excluded_keys=['IO'])
        key_reads = 'PE' if params.read_type == 'illumina' else 'SE'
        quast_.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), read_type=key_reads))
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
        FASTA = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA
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
        FASTA =  Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_FASTA,
        DIR = rules.quast_quast.output.DIR,
        INFORMS_quast = rules.quast_quast.output.INFORMS,
        INFORMS_assembler = Path(config['working_dir']) / assembly_spades.OUTPUT_ASSEMBLY_INFORMS,
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
        ]

        # Parse QUAST report
        path_tsv = SnakemakeUtils.load_object(Path(input.TSV))[0].path
        data_quast = {}
        with path_tsv.open() as handle:
            for line in handle.readlines():
                print(line.strip())
                key, value = line.strip().split('\t')
                data_quast[key] = value

        # Create TSV output
        with open(output.TSV, 'w') as handle:
            for row in keys_kept:
                handle.write(f"assembly_{row['name']}\t{data_quast[row['key']]}")
                handle.write('\n')
