import os

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import SNAKEFILE_SEQUENCE_TYPING_BLAST, SNAKEFILE_SEQUENCE_TYPING_SRST2
from camel.resources.snakefile.sequence_typing import OUTPUT_TYPING_REPORT, OUTPUT_TYPING_REPORT_EMPTY, \
    OUTPUT_TYPING_SUMMARY, OUTPUT_TYPING_HITS

##################
#  Configuration #
##################
camel = Camel.get_instance()
SCHEMES = config['sequence_typing']
SCHEME_METADATA = {name: SequenceTypingUtils.parse_scheme_metadata(path) for name, path in SCHEMES.items()}
loci_by_scheme_by_type = {
    name: {
        'DNA': [locus['name_valid'] for locus in metadata['loci'] if locus['type'] == 'DNA'],
        'peptide': [locus['name_valid'] for locus in metadata['loci'] if locus['type'] == 'peptide']
    } for name, metadata in SCHEME_METADATA.items()
}

##############################
# Allele detection workflows #
##############################
include: SNAKEFILE_SEQUENCE_TYPING_BLAST
include: SNAKEFILE_SEQUENCE_TYPING_SRST2

#########
# Rules #
#########
rule Typing_extract_schema_info:
    """
    Extracts the metadata for a scheme.
    """
    input:
        lambda wildcards: SCHEMES[wildcards.scheme]
    output:
        INFORMS=os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-locus_set.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'typing', '{scheme}')
    run:
        from camel.app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
        locus_set_manager = LocusSetManager(camel)
        locus_set_manager.add_input_files({'DIR': [ToolIODirectory(input[0])]})
        step = Step(rule, locus_set_manager, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_object(locus_set_manager.informs, output.INFORMS)

rule Typing_pickle_dump_sequence_type_definitions:
    """
    Retrieves the sequence type definitions and converts them to CAMEL IO pickles.
    """
    input:
        lambda wildcards: SCHEMES[wildcards.scheme]
    output:
        os.path.join(config['working_dir'], 'typing', '{scheme}', 'tsv-profiles.io')
    run:
        SnakemakeUtils.dump_object([ToolIOFile(os.path.join(input[0], 'profiles.tsv'))], output[0])

rule Typing_get_hits:
    """
    Selects the hits output based on the detection method in the config.
    """
    input:
        hits_nucl=lambda wildcards: os.path.join(config['working_dir'], OUTPUT_TYPING_HITS.format(scheme=wildcards.scheme, locus_type='DNA', detection_method=config['detection_method'])) if (len(loci_by_scheme_by_type[wildcards.scheme]['DNA']) > 0) else [],
        hits_pept=lambda wildcards: os.path.join(config['working_dir'], OUTPUT_TYPING_HITS.format(scheme=wildcards.scheme, locus_type='peptide', detection_method='blast')) if (len(loci_by_scheme_by_type[wildcards.scheme]['peptide']) > 0) else []
    output:
        hits_nucl=os.path.join(config['working_dir'], 'typing', '{scheme}', 'DNA', 'hits.io'),
        hits_pept=os.path.join(config['working_dir'], 'typing', '{scheme}', 'peptide', 'hits.io')
    run:
        import shutil
        for key in ('hits_nucl', 'hits_pept'):
            data = input[key]
            if len(data) > 0:
                shutil.copyfile(data, output.get(key))
            else:
                SnakemakeUtils.dump_object([], output.get(key))

rule Typing_export_hits_tabular:
    """
    Creates a tabular output for the detected hits.
    """
    input:
        hits=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', 'hits.io')
    output:
        TSV=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', 'tabular', 'tsv.io')
    params:
        working_dir=os.path.join(config['working_dir'], 'typing', '{scheme}', '{locus_type}', 'tabular'),
        sample_name=FileSystemHelper.make_valid(config['sample_name']),
        scheme=lambda wildcards: FileSystemHelper.make_valid(wildcards.scheme),
        locus_type=lambda wildcards: wildcards.locus_type
    run:
        hits = SnakemakeUtils.load_object(input.hits)
        output_file = os.path.join(
            params.working_dir, f'typing-{params.scheme}-{params.locus_type}-{params.sample_name}.tsv')
        if len(hits) == 0:
            SnakemakeUtils.dump_object([], output.TSV)
        else:
            with open(output_file, 'w') as handle_out:
                handle_out.write('\t'.join(hits[0].value.get_table_column_names()))
                handle_out.write('\n')
                for h in hits:
                    handle_out.write(h.value.to_table_row())
                    handle_out.write('\n')
            SnakemakeUtils.dump_object([ToolIOFile(output_file)], output.TSV)

rule Typing_detect_sequence_type:
    """
    Detects the sequence type based on the detected alleles.
    """
    input:
         hits_nucl=os.path.join(config['working_dir'], 'typing', '{scheme}', 'DNA', 'hits.io'),
         hits_pept=os.path.join(config['working_dir'], 'typing', '{scheme}', 'peptide', 'hits.io'),
         TSV=os.path.join(config['working_dir'], 'typing', '{scheme}', 'tsv-profiles.io')
    output:
        INFORMS=os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-st.io')
    params:
        running_dir=os.path.join(config['working_dir'], 'typing', '{scheme}')
    run:
        from camel.app.tools.pipelines.sequence_typing.sequencetypedetector import SequenceTypeDetector
        sequence_type_detector = SequenceTypeDetector(camel)
        SnakemakeUtils.add_pickle_inputs(sequence_type_detector, input)
        step = Step(rule, sequence_type_detector, camel, params.running_dir, config)
        sequence_type_detector.update_parameters(allele_wildcard='N', allele_absent_symbol='0')
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(sequence_type_detector, output)

rule Typing_add_allele_page_url:
    """
    This steps add the locus url to the detected hits.
    """
    input:
        hits_nucl=os.path.join(config['working_dir'], 'typing', '{scheme}', 'DNA', 'hits.io'),
        hits_pept=os.path.join(config['working_dir'], 'typing', '{scheme}', 'peptide', 'hits.io'),
        INFORMS_scheme = os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-locus_set.io')
    output:
        hits_nucl=os.path.join(config['working_dir'], 'typing', '{scheme}', 'DNA', 'hits-url.io'),
        hits_pept=os.path.join(config['working_dir'], 'typing', '{scheme}', 'peptide', 'hits-url.io')
    run:
        for key in ('nucl', 'pept'):
            hits = SnakemakeUtils.load_object(input.get(f'hits_{key}'))

            # Load the informs
            metadata_by_locus_name = SnakemakeUtils.load_object(input.INFORMS_scheme)['loci'].metadata_by_locus_name

            # Add the allele url to the hit
            for hit in hits:
                locus_key = FileSystemHelper.make_valid(hit.value.locus).lower()
                locus_metadata = metadata_by_locus_name[key]
                hit.value.set_allele_page_url_template(locus_metadata.get('allele_page_url'))

            # Export hits
            SnakemakeUtils.dump_object(hits, output.get(f'hits_{key}'))

rule Typing_create_report:
    """
    Creates a report with the sequence typing output.
    """
    input:
        TSV_nucl=os.path.join(config['working_dir'], 'typing', '{scheme}', 'DNA', 'tabular', 'tsv.io'),
        TSV_pept=os.path.join(config['working_dir'], 'typing', '{scheme}', 'peptide', 'tabular', 'tsv.io'),
        INFORMS_scheme=os.path.join(config['working_dir'], 'typing' , '{scheme}', 'informs-locus_set.io'),
        INFORMS_ST=lambda wildcards: os.path.join(config['working_dir'], 'typing',  wildcards.scheme, 'informs-st.io') if SequenceTypingUtils.has_profiles(SCHEMES, wildcards.scheme) else [],
        hits_nucl=os.path.join(config['working_dir'], 'typing', '{scheme}', 'DNA', 'hits-url.io'),
        hits_pept=os.path.join(config['working_dir'], 'typing', '{scheme}', 'peptide', 'hits-url.io')
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_TYPING_REPORT)
    params:
        running_dir=os.path.join(config['working_dir'], 'typing', '{scheme}'),
        sample_name=config['sample_name']
    run:
        from camel.app.tools.pipelines.sequence_typing.htmlreportertyping import HtmlReporterTyping
        reporter = HtmlReporterTyping(camel)
        if len(input.INFORMS_ST) != 0:
           reporter.add_input_informs({'ST': SnakemakeUtils.load_object(input.INFORMS_ST)})
        SnakemakeUtils.add_pickle_inputs(reporter, input, excluded_keys=['INFORMS_ST'])
        reporter.add_input_files({'VAL_SAMPLE': [ToolIOValue(params.sample_name)]})
        step = Step(rule, reporter, camel, params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule Typing_create_report_empty:
    """
    Creates an empty sequence typing report when the analysis is disabled.
    """
    input:
        INFORMS_Scheme=os.path.join(config['working_dir'], 'typing', '{scheme}', 'informs-locus_set.io')
    output:
        VAL_HTML=os.path.join(config['working_dir'], OUTPUT_TYPING_REPORT_EMPTY)
    run:
        informs = SnakemakeUtils.load_object(input.INFORMS_Scheme)
        section = HtmlReportSection(informs['title'], 3)
        section.add_paragraph('Analysis disabled')
        SnakemakeUtils.dump_object([ToolIOValue(section)], output[0])

rule Typing_dump_summary_info:
    """
    Dumps the summary information in tabular format.
    """
    input:
        hits_nucl=os.path.join(config['working_dir'], 'typing', '{scheme}', 'DNA', 'hits-url.io'),
        hits_pept=os.path.join(config['working_dir'], 'typing', '{scheme}', 'peptide', 'hits-url.io'),
        INFORMS_ST=lambda wildcards: os.path.join(config['working_dir'], 'typing', wildcards.scheme, 'informs-st.io') if SequenceTypingUtils.has_profiles(SCHEMES, wildcards.scheme) else []
    output:
        os.path.join(config['working_dir'], OUTPUT_TYPING_SUMMARY)
    params:
        scheme_name=lambda wildcards: wildcards.scheme
    run:
        if len(input.INFORMS_ST) == 0:
            st_metadata = []
        else:
            st_metadata = SnakemakeUtils.load_object(input.INFORMS_ST)['sequence_type'].metadata
        hits = SnakemakeUtils.load_object(input.hits_nucl) + SnakemakeUtils.load_object(input.hits_pept)
        with open(output[0], 'w') as handle:
            for k, v in st_metadata:
                handle.write(f'{params.scheme_name}-{k}\t{v}')
                handle.write('\n')
            for hit in hits:
                key = '{}-{}'.format(params.scheme_name, hit.value.locus)
                allele_id = hit.value.to_table_row(separator=',')
                handle.write(f'{key}\t{allele_id}')
                handle.write('\n')
