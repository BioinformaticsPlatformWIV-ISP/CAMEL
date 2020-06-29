from pathlib import Path

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import sequence_typing
from camel.resources.snakefile.sequence_typing import OUTPUT_TYPING_INFORMS


##################
#  Configuration #
##################
camel = Camel.get_instance()
SCHEME_DATA = config['sequence_typing']
SCHEME_METADATA = {name: SequenceTypingUtils.parse_scheme_metadata(data['path']) for name, data in SCHEME_DATA.items()}
loci_by_scheme_by_type = {
    name: {
        'DNA': [locus['name_valid'] for locus in metadata['loci'] if locus['type'] == 'DNA'],
        'peptide': [locus['name_valid'] for locus in metadata['loci'] if locus['type'] == 'peptide']
    } for name, metadata in SCHEME_METADATA.items()
}

##############################
# Allele detection workflows #
##############################
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING_BLAST
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING_SRST2
include: sequence_typing.SNAKEFILE_SEQUENCE_TYPING_KMA

#########
# Rules #
#########
rule typing_extract_schema_info:
    """
    Extracts the metadata for a scheme.
    """
    input:
        lambda wildcards: SCHEME_DATA[wildcards.scheme]['path']
    output:
        INFORMS = Path(config['working_dir']) / 'typing' / '{scheme}' /'informs-locus_set.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'typing' / wildcards.scheme
    run:
        from camel.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
        locus_set_manager = LocusSetManager(camel)
        locus_set_manager.add_input_files({'DIR': [ToolIODirectory(input[0])]})
        step = Step(rule, locus_set_manager, camel, params.running_dir, config, wildcards)
        step.run_step()
        SnakemakeUtils.dump_object(locus_set_manager.informs, output.INFORMS)

rule typing_pickle_profiles:
    """
    Retrieves the sequence type definitions and converts them to CAMEL IO pickles.
    """
    input:
        DIR = lambda wildcards: SCHEME_DATA[wildcards.scheme]['path']
    output:
        TSV = Path(config['working_dir']) / 'typing' / '{scheme}' / 'tsv-profiles.io'
    run:
        from camel.app.io.tooliofile import ToolIOFile
        SnakemakeUtils.dump_object([ToolIOFile(Path(input.DIR) / 'profiles.tsv')], output.TSV)

rule typing_get_hits:
    """
    Selects the hits output based on the detection method in the config.
    """
    input:
        HITS_NUCL = lambda wildcards: str(Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_HITS).format(
            scheme=wildcards.scheme,
            locus_type='DNA',
            detection_method=SequenceTypingUtils.get_detection_method(config, wildcards.scheme))) if (
                len(loci_by_scheme_by_type[wildcards.scheme]['DNA']) > 0) else [],
        HITS_PEPT = lambda wildcards: str(Path(config['working_dir']) / str(sequence_typing.OUTPUT_TYPING_HITS).format(
            scheme=wildcards.scheme,
            locus_type='peptide',
            detection_method='blast')) if (len(loci_by_scheme_by_type[wildcards.scheme]['peptide']) > 0) else []
    output:
        HITS_NUCL = Path(config['working_dir']) / 'typing' / '{scheme}' / 'DNA' / 'hits.io',
        HITS_PEPT = Path(config['working_dir']) / 'typing' / '{scheme}' / 'peptide' / 'hits.io'
    run:
        import shutil
        for key in ('HITS_NUCL', 'HITS_PEPT'):
            data = input[key]
            if len(data) > 0:
                shutil.copyfile(data, output.get(key))
            else:
                SnakemakeUtils.dump_object([], output.get(key))

rule typing_dump_commands:
    input:
        IO = lambda wildcards: [str(Path(config['working_dir']) / 'typing' / wildcards.scheme / type_ / loci_by_scheme_by_type[wildcards.scheme][type_][0] / f'informs-{SequenceTypingUtils.get_detection_method(config, wildcards.scheme)}.io') for type_ in (
            'DNA', 'peptide') if len(loci_by_scheme_by_type[wildcards.scheme][type_]) > 0],
        INFORMS_scheme = Path(config['working_dir']) / 'typing' / '{scheme}' /'informs-locus_set.io'
    output:
        INFORMS = Path(config['working_dir']) / OUTPUT_TYPING_INFORMS
    params:
        scheme = lambda wildcards: wildcards.scheme
    run:
        informs_scheme = SnakemakeUtils.load_object(input.INFORMS_scheme)
        informs_commands = [SnakemakeUtils.load_object(io) for io in input.IO]
        for informs in informs_commands:
            informs['_tag'] = f"Typing - {informs_scheme['title']}"
        SnakemakeUtils.dump_object(informs_commands, output.INFORMS)

rule typing_export_hits_tabular:
    """
    Creates a tabular output for the detected hits.
    """
    input:
        hits = Path(config['working_dir']) / 'typing' / '{scheme}' / '{locus_type}' / 'hits.io'
    output:
        TSV = Path(config['working_dir']) / 'typing' / '{scheme}' / '{locus_type}' / 'tabular' / 'tsv.io'
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / 'typing' / wildcards.scheme / wildcards.locus_type / 'tabular',
        sample_name = FileSystemHelper.make_valid(config['sample_name']),
        scheme = lambda wildcards: FileSystemHelper.make_valid(wildcards.scheme),
        locus_type = lambda wildcards: wildcards.locus_type
    run:
        hits = SnakemakeUtils.load_object(input.hits)
        output_file = params.working_dir / f'typing-{params.scheme}-{params.locus_type}-{params.sample_name}.tsv'
        if len(hits) == 0:
            SnakemakeUtils.dump_object([], output.TSV)
        else:
            with output_file.open('w') as handle_out:
                handle_out.write('\t'.join(hits[0].value.table_column_names()))
                handle_out.write('\n')
                for h in hits:
                    handle_out.write('\t'.join(h.value.to_table_row()))
                    handle_out.write('\n')
            SnakemakeUtils.dump_object([ToolIOFile(output_file)], output.TSV)

rule typing_detect_sequence_type:
    """
    Detects the sequence type based on the detected alleles.
    """
    input:
         hits_nucl = rules.typing_get_hits.output.HITS_NUCL,
         hits_pept = rules.typing_get_hits.output.HITS_PEPT,
         TSV = rules.typing_pickle_profiles.output.TSV
    output:
        INFORMS = Path(config['working_dir']) / 'typing' / '{scheme}' / 'informs-st.io'
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'typing' / wildcards.scheme
    run:
        from camel.app.tools.pipelines.sequence_typing.sequencetypedetector import SequenceTypeDetector
        sequence_type_detector = SequenceTypeDetector(camel)
        SnakemakeUtils.add_pickle_inputs(sequence_type_detector, input)
        step = Step(rule, sequence_type_detector, camel, params.running_dir, config, wildcards)
        sequence_type_detector.update_parameters(allele_wildcard='N', allele_absent_symbol='0')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sequence_type_detector, output)

rule typing_get_cgmlst_stats:
    """
    Retrieves the cgMLST stats for the given scheme.
    Only DNA loci are considered.
    """
    input:
        HITS = rules.typing_get_hits.output.HITS_NUCL
    output:
        INFORMS = Path(config['working_dir']) / 'typing' / '{scheme}' / 'stats' / 'informs.io'
    params:
        scheme_name = lambda wildcards: wildcards.scheme
    run:
        all_hits = SnakemakeUtils.load_object(input.HITS)
        nb_perfect = len([v for v in all_hits if v.value.is_perfect_hit()])
        SnakemakeUtils.dump_object(
            {'hits_found': nb_perfect,
             'nb_of_loci': len(all_hits),
             'scheme_name': params.scheme_name},
            output.INFORMS)

rule typing_add_allele_page_url:
    """
    This steps add the locus url to the detected hits.
    """
    input:
        HITS_NUCL = rules.typing_get_hits.output.HITS_NUCL,
        HITS_PEPT = rules.typing_get_hits.output.HITS_PEPT,
        INFORMS_scheme = rules.typing_extract_schema_info.output.INFORMS
    output:
        HITS_NUCL = Path(config['working_dir']) / 'typing' / '{scheme}' / 'DNA' / 'hits-url.io',
        HITS_PEPT = Path(config['working_dir']) / 'typing' / '{scheme}' / 'peptide' / 'hits-url.io'
    run:
        for key in ('NUCL', 'PEPT'):
            hits = SnakemakeUtils.load_object(input.get(f'HITS_{key}'))

            # Load the informs
            metadata_by_locus_name = SnakemakeUtils.load_object(input.INFORMS_scheme)['loci'].metadata_by_locus_name

            # Add the allele url to the hit
            for hit in hits:
                locus_key = FileSystemHelper.make_valid(hit.value.locus).lower()
                locus_metadata = metadata_by_locus_name[locus_key]
                hit.value.set_allele_page_url_template(locus_metadata.get('allele_page_url'))

            # Export hits
            SnakemakeUtils.dump_object(hits, output.get(f'HITS_{key}'))

rule typing_create_report:
    """
    Creates a report with the sequence typing output.
    """
    input:
        TSV_nucl = rules.typing_export_hits_tabular.output.TSV.format(locus_type='DNA', scheme='{scheme}'),
        TSV_pept = rules.typing_export_hits_tabular.output.TSV.format(locus_type='peptide', scheme='{scheme}'),
        INFORMS_scheme = rules.typing_extract_schema_info.output.INFORMS,
        INFORMS_ST = lambda wildcards: rules.typing_detect_sequence_type.output.INFORMS if SequenceTypingUtils.has_profiles(SCHEME_DATA, wildcards.scheme) else [],
        hits_nucl = rules.typing_get_hits.output.HITS_NUCL,
        hits_pept = rules.typing_get_hits.output.HITS_PEPT
    output:
        VAL_HTML = Path(config['working_dir']) / sequence_typing.OUTPUT_TYPING_REPORT
    params:
        running_dir = lambda wildcards: Path(config['working_dir']) / 'typing' / wildcards.scheme,
        sample_name = config['sample_name'],
        detection_method = lambda wildcards: SequenceTypingUtils.get_detection_method(config, wildcards.scheme)
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.sequence_typing.htmlreportertyping import HtmlReporterTyping
        reporter = HtmlReporterTyping(camel)
        if len(input.INFORMS_ST) != 0:
           reporter.add_input_informs({'ST': SnakemakeUtils.load_object(input.INFORMS_ST)})
        SnakemakeUtils.add_pickle_inputs(reporter, input, excluded_keys=['INFORMS_ST'])
        reporter.add_input_files({'VAL_SAMPLE': [ToolIOValue(params.sample_name)]})
        if params.detection_method != config['detection_method']:
            reporter.update_parameters(forced_detection_method=str(params.detection_method))
        step = Step(rule, reporter, camel, params.running_dir, config, wildcards)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule typing_create_report_empty:
    """
    Creates an empty sequence typing report when the analysis is disabled.
    """
    input:
        INFORMS_Scheme = rules.typing_extract_schema_info.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / sequence_typing.OUTPUT_TYPING_REPORT_EMPTY
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        informs = SnakemakeUtils.load_object(input.INFORMS_Scheme)
        section = HtmlReportSection(informs['title'], 3)
        section.add_paragraph('Analysis disabled')
        SnakemakeUtils.dump_object([ToolIOValue(section)], output[0])

rule typing_dump_summary_info:
    """
    Dumps the summary information in tabular format.
    """
    input:
        HITS_NUCL = rules.typing_add_allele_page_url.output.HITS_NUCL,
        HITS_PEPT = rules.typing_add_allele_page_url.output.HITS_PEPT,
        INFORMS_ST = lambda wildcards: rules.typing_detect_sequence_type.output.INFORMS if SequenceTypingUtils.has_profiles(SCHEME_DATA, wildcards.scheme) else []
    output:
        TSV = Path(config['working_dir']) / sequence_typing.OUTPUT_TYPING_SUMMARY
    params:
        scheme_name = lambda wildcards: wildcards.scheme
    run:
        if len(input.INFORMS_ST) == 0:
            st_metadata = []
        else:
            st_metadata = SnakemakeUtils.load_object(input.INFORMS_ST)['sequence_type'].metadata
        hits = SnakemakeUtils.load_object(input.HITS_NUCL) + SnakemakeUtils.load_object(input.HITS_PEPT)
        with open(output.TSV, 'w') as handle:
            for k, v in st_metadata:
                handle.write(f'{params.scheme_name}-{k}\t{v}')
                handle.write('\n')
            for hit in hits:
                key = '{}-{}'.format(params.scheme_name, hit.value.locus)
                allele_id = ','.join(hit.value.to_table_row())
                handle.write(f'{key}\t{allele_id}')
                handle.write('\n')
