from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import medaka_polishing


rule medaka_polishing_map_ont_reads:
    """
    Maps the ONT reads to the Flye assembly.
    """
    input:
        FQ = Path(config['working_dir']) / 'fq_dict.io',
        FASTA = lambda wildcards: Path(config['working_dir']) / str(medaka_polishing.INPUT_ASSEMBLY_FASTA).format(assembly_type=wildcards.assembly_type)
    output:
        BAM = Path(config['working_dir']) / 'medaka' / '{assembly_type}' / 'minimap2' / 'bam.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'medaka' / wildcards.assembly_type / 'minimap2'
    threads: 8
    run:
        from camel.app.components.pipelines import pipeutils
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex

        # Minimap2
        minimap2 = Minimap2Mapping(Camel.get_instance())
        SnakemakeUtils.add_pickle_input(minimap2, 'FASTA', Path(str(input.FASTA)))
        minimap2.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.FQ), key_se='FASTQ', read_type='SE'))

        # Initialize tools
        samtools_view = SamtoolsView(Camel.get_instance())
        samtools_sort = SamtoolsSort(Camel.get_instance())
        samtools_sort.update_parameters(threads=threads)
        pipeutils.run_as_pipe([minimap2, samtools_view, samtools_sort], Path(str(params.dir_)))

        # Create BAM index
        samtools_index = SamtoolsIndex(Camel.get_instance())
        samtools_index.update_parameters(generate_bai_index=True)
        samtools_index.add_input_files({'BAM': samtools_sort.tool_outputs['BAM']})
        samtools_index.run(Path(str(params.dir_)))

        # Export output
        SnakemakeUtils.dump_tool_output(samtools_index, 'BAM', Path(output.BAM))

rule medaka_polishing_medaka_consensus:
    """
    Runs the medaka consensus algorithm.
    """
    input:
        BAM = rules.medaka_polishing_map_ont_reads.output.BAM
    output:
        HDF = Path(config['working_dir']) / 'medaka' / '{assembly_type}' / 'consensus' / 'raw_hdf.io',
        INFORMS = Path(config['working_dir']) / 'medaka' / '{assembly_type}' / 'consensus' / 'commands-consensus.io'
    params:
        dir_ =  lambda wildcards: Path(config['working_dir']) / 'medaka' / wildcards.assembly_type / 'consensus',
        medaka_options = config.get('polishing', {}).get('medaka', {}).get('consensus', {})
    threads: 8
    run:
        from camel.app.tools.medaka.medakaconsensus import MedakaConsensus
        medaka = MedakaConsensus(Camel.get_instance())
        SnakemakeUtils.add_pickle_input(medaka, 'BAM', Path(input.BAM))
        medaka.update_parameters(**params.medaka_options)
        medaka.update_parameters(threads=threads)
        step = Step(str(rule), medaka, Camel.get_instance(), Path(str(params.dir_)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(medaka, output)

rule medaka_polishing_medaka_stitch:
    """
    Runs the medaka stitch algorithm.
    """
    input:
        HDF = rules.medaka_polishing_medaka_consensus.output.HDF,
        FASTA = Path(config['working_dir']) / medaka_polishing.INPUT_ASSEMBLY_FASTA
    output:
        FASTA = Path(config['working_dir']) / 'medaka' / '{assembly_type}' / 'stitch' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / 'medaka' / '{assembly_type}' / 'stitch' / 'commands-stitch.io'
    params:
        dir_ = lambda wildcards: Path(config['working_dir']) / 'medaka' / wildcards.assembly_type / 'stitch',
        medaka_options = config.get('polishing', {}).get('medaka', {}).get('stitch', {})
    threads: 8
    run:
        from camel.app.tools.medaka.medakastitch import MedakaStitch
        medaka = MedakaStitch(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(medaka, input)
        medaka.update_parameters(**params.medaka_options, threads=threads)
        step = Step(str(rule), medaka, Camel.get_instance(), Path(str(params.dir_)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(medaka, output)

rule medaka_polishing_empty_report:
    """
    Creates an empty report for the gene detection when plasmidSPAdes assembly fails.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / medaka_polishing.OUTPUT_ASSEMBLY_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Medaka polishing', Path(output.VAL_HTML))
