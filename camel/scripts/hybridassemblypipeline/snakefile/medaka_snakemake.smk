import pickle
from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step

camel = Camel.get_instance()

rule map_reads_to_assembly:
    input:
        FASTQ = Path(config['working_dir']) / 'trimming' / 'ont' / '{}_SE.fastq.gz'.format(config['name']),
        FASTA = Path(config['working_dir']) / 'assembly_flye' / 'filtering' / 'assembly_filtered.fasta'
    output:
        SAM =  Path(config['working_dir']) / 'mapping' / 'reads_to_flye_assembly.sam'
    params:
        working_dir = Path(config['working_dir'] / 'mapping')
    threads: 4
    run:
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        minimap = Minimap2Mapping(camel)
        minimap.add_input_files({'FASTQ': [ToolIOFile(Path(input.FASTQ))], 'FASTA': [ToolIOFile(Path(input.FASTA))]})
        minimap.update_parameters(output_filename='reads_to_flye_assembly.sam')
        step = Step(str(rule),minimap,camel,params.working_dir,config)
        step.run_step()

rule sam_to_bam:
    input:
        SAM = rules.map_reads_to_assembly.output.SAM
    output:
        BAM = Path(config['working_dir']) / 'mapping' / 'reads_to_flye_assembly.bam'
    params:
        working_dir = Path(config['working_dir'] / 'mapping')
    threads: 4
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        samtools = SamtoolsView(camel)
        samtools.add_input_files({'SAM': [ToolIOFile(Path(input.SAM))]})
        samtools.update_parameters(output_filename='reads_to_flye_assembly.bam')
        step = Step(str(rule),samtools,camel,params.working_dir,config)
        step.run_step()

rule sort_bam:
    input:
        BAM = rules.sam_to_bam.output.BAM
    output:
        BAM = Path(config['working_dir']) / 'mapping' / 'reads_to_flye_assembly.sorted.bam'
    params:
        working_dir = Path(config['working_dir'] / 'mapping')
    threads: 4
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools = SamtoolsSort(camel)
        samtools.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        samtools.update_parameters(output_filename='reads_to_flye_assembly.sorted.bam')
        step = Step(str(rule),samtools,camel,params.working_dir,config)
        step.run_step()

rule index_bam:
    input:
        BAM = rules.sort_bam.output.BAM
    output:
        BAI = Path(config['working_dir']) / 'mapping' / 'reads_to_flye_assembly.sorted.bam.bai'
    params:
        working_dir = Path(config['working_dir'] / 'mapping')
    threads: 4
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        samtools = SamtoolsIndex(camel)
        samtools.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        samtools.update_parameters(generate_bai_index=True)
        step = Step(str(rule),samtools,camel,params.working_dir,config)
        step.run_step()

rule medaka_consensus:
    """
    Runs the medaka consensus algorithm.
    """
    input:
        BAM = rules.sort_bam.output.BAM,
        BAI = rules.index_bam.output.BAI
    output:
        HDF = Path(config['working_dir']) / 'medaka' / 'raw.hdf',
        INFORMS = Path(config['working_dir']) / 'medaka' / 'commands-consensus.io'
    params:
        working_dir= Path(config['working_dir']) / 'medaka',
        medaka_options = config.get('polishing', {}).get('medaka', {}).get('consensus', {})
    threads: 8
    run:
        from camel.app.tools.medaka.medakaconsensus import MedakaConsensus
        medaka = MedakaConsensus(camel)
        medaka.add_input_files({'BAM': [ToolIOFile(Path(input.BAM))]})
        medaka.update_parameters(**params.medaka_options)
        medaka.update_parameters(threads=threads)
        step = Step(str(rule), medaka, camel, params.working_dir, config)
        step.run_step()
        with open(output.INFORMS, 'wb') as handle:
            pickle.dump(medaka.informs, handle)

rule medaka_stitch:
    """
    Runs the medaka stitch algorithm.
    """
    input:
        HDF = rules.medaka_consensus.output.HDF,
        FASTA = Path(config['working_dir']) / 'assembly_flye' / 'filtering' / 'assembly_filtered.fasta'
    output:
        FASTA = Path(config['working_dir']) / 'medaka' / 'consensus.fasta',
        INFORMS = Path(config['working_dir']) / 'medaka' / 'commands-stitch.io'
    params:
        working_dir= Path(config['working_dir']) / 'medaka',
        medaka_options = config.get('polishing',{}).get('medaka',{}).get('stitch',{})
    threads: 8
    run:
        from camel.app.tools.medaka.medakastitch import MedakaStitch
        medaka = MedakaStitch(camel)
        medaka.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA))], 'HDF': [ToolIOFile(Path(input.HDF))]})
        medaka.update_parameters(**params.medaka_options)
        medaka.update_parameters(threads=threads)
        step = Step(str(rule), medaka, camel, params.working_dir, config)
        step.run_step()
        with open(output.INFORMS, 'wb') as handle:
            pickle.dump(medaka.informs, handle)
