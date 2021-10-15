from pathlib import Path

from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import variant_calling, gene_detection
from camel.scripts.shigellapipeline.snakefile import subspecies_identification

rule subspecies_identification_speG_depth:
    """
    Calculates the average depth of coverage of the speG gene.
    """
    input:
        BAM = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_BAM
    output:
        INFORMS = Path(config['working_dir']) / 'subspecies_identification' / 'informs-depth.io'
    params:
        bed_file = config['subspecies_identification']['bed_speG'],
        running_dir = Path(config['working_dir']) / 'subspecies_identification'
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth(camel)
        SnakemakeUtils.add_pickle_inputs(samtools_depth, input)
        samtools_depth.add_input_files({'BED': [ToolIOFile(Path(params.bed_file))]})
        step = Step(rule, samtools_depth, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(samtools_depth, output)

rule subspecies_identification_detect_species:
    """
    Detects the species (E. coli / Shigella).
    """
    input:
        VAL_hits = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_ALL_HITS).format(db='subspecies_identification'),
        INFORMS_speG_depth = rules.subspecies_identification_speG_depth.output.INFORMS,
        VCF = Path(config['working_dir']) / variant_calling.OUTPUT_VARIANT_CALLING_UNFILTERED_VCF
    output:
        INFORMS = Path(config['working_dir']) / 'subspecies_identification' / 'species' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'subspecies_identification' / 'species'
    run:
        from camel.app.tools.pipelines.shigella.speciesdetector import SpeciesDetector
        species_detector = SpeciesDetector(camel)
        SnakemakeUtils.add_pickle_inputs(species_detector, input)
        step = Step(rule, species_detector, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(species_detector, output)

rule subspecies_identification_report_species:
    """
    Creates a report with the species identification.
    """
    input:
        VAL_hits = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_ALL_HITS).format(db='subspecies_identification'),
        INFORMS_columns=Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_COLUMNS).format(db='subspecies_identification'),
        INFORMS_species = rules.subspecies_identification_detect_species.output.INFORMS
    output:
        VAL_HTML_species = Path(config['working_dir']) / 'subspecies_identification' / 'report' / 'html-species.io'
    params:
        running_dir = Path(config['working_dir']) / 'subspecies_identification' / 'report'
    run:
        from camel.app.tools.pipelines.shigella.speciesreporter import SpeciesReporter
        reporter = SpeciesReporter(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule subspecies_identification_report_species_empty:
    """
    Creates an empty report when (sub)species identification is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / subspecies_identification.OUTPUT_SPECIES_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Species identification', Path(output.VAL_HTML))

rule subspecies_identification_detect_subspecies:
    """
    Detects the Shigella subspecies.
    """
    input:
        VAL_hits = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_ALL_HITS).format(db='subspecies_identification'),
    output:
        INFORMS = Path(config['working_dir']) / subspecies_identification.OUTPUT_SPECIES_SUBSPECIES_INFORMS
    params:
        running_dir = Path(config['working_dir']) / 'subspecies_identification' / 'subspecies'
    run:
        from camel.app.tools.pipelines.shigella.subspeciesdetector import SubspeciesDetector
        subspecies_detector = SubspeciesDetector(camel)
        SnakemakeUtils.add_pickle_inputs(subspecies_detector, input)
        step = Step(rule, subspecies_detector, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(subspecies_detector, output)

rule subspecies_identification_report_subspecies:
    """
    Creates a report with the subspecies identification.
    """
    input:
        VAL_hits = Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_ALL_HITS).format(db='subspecies_identification'),
        INFORMS_columns=Path(config['working_dir']) / str(gene_detection.OUTPUT_GENE_DETECTION_COLUMNS).format(db='subspecies_identification'),
        INFORMS_subspecies = rules.subspecies_identification_detect_subspecies.output.INFORMS
    output:
        VAL_HTML_subspecies = Path(config['working_dir']) / 'subspecies_identification' / 'report' / 'html-subspecies.io'
    params:
        running_dir = Path(config['working_dir']) / 'subspecies_identification' / 'report'
    run:
        from camel.app.tools.pipelines.shigella.subspeciesreporter import SubspeciesReporter
        reporter = SubspeciesReporter(camel)
        SnakemakeUtils.add_pickle_inputs(reporter, input)
        step = Step(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule subspecies_identification_report_subspecies_empty:
    """
    Creates an empty report when (sub)species identification is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / subspecies_identification.OUTPUT_SUBSPECIES_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Subspecies identification', Path(output.VAL_HTML))

rule subspecies_identification_dump_summary_info:
    """
    Dumps the summary information for the subspecies identification assay.
    """
    input:
        INFORMS_subspecies = rules.subspecies_identification_detect_subspecies.output.INFORMS,
        INFORMS_species = rules.subspecies_identification_detect_species.output.INFORMS
    output:
        TSV = Path(config['working_dir']) / subspecies_identification.OUTPUT_SUBSPECIES_SUMMARY
    run:
        informs_subspecies = SnakemakeUtils.load_object(Path(input.INFORMS_subspecies))
        informs_species = SnakemakeUtils.load_object(Path(input.INFORMS_species))
        output_data = [
            ['detected_species', informs_species['detected_species']],
            ['speG_present', informs_species['speG_present']],
            ['speG_indel_present', informs_species['speG_indel_present']],
            ['ipaH_present', informs_species['ipaH_present']],
            ['detected_subspecies', informs_subspecies['detected_subspecies']],
            ['detected_loci_subspecies', informs_subspecies['detected_loci_subspecies']],
        ]
        with open(output.TSV, 'w') as handle:
            for key, value in output_data:
                handle.write(f'{key}\t{value}\n')

rule subspecies_identification_export_database_info:
    """
    Creates the report section with the database sequences used in the subspecies identification.
    """
    input:
        FASTA_subspecies = Path(config['working_dir']) / 'gene_detection' / 'subspecies_identification' / 'db_manager' / 'fasta.io',
        FASTA_flexneri = Path(config['working_dir']) / 'gene_detection' / 'flexneri_type' / 'db_manager' / 'fasta.io'
    output:
        Path(config['working_dir']) / 'subspecies_identification' / 'report' / 'html-db.io'
    params:
        fasta_promotor = config['flexneri_type']['fasta_gtr_promotor']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.components.html.htmlreportsection import HtmlReportSection
        section = HtmlReportSection('Database information')
        relative_path_subspecies = section.add_file(
            SnakemakeUtils.load_object(Path(input.FASTA_subspecies))[0].path, Path('subspecies') / 'subspecies.fasta')
        section.add_link_to_file('Subspecies sequences (FASTA)', relative_path_subspecies)

        relative_path_flexneri = section.add_file(
            SnakemakeUtils.load_object(Path(input.FASTA_subspecies))[0].path, Path('subspecies', 'flexneri.fasta'))
        section.add_link_to_file('Flexneri type sequences (FASTA)', relative_path_flexneri)

        relative_path_gtr = section.add_file(
            Path(params.fasta_promotor), Path('subspecies', 'gtr_promotor.fasta'))
        section.add_link_to_file('gtr promotor (FASTA)', relative_path_gtr)

        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output[0]))
