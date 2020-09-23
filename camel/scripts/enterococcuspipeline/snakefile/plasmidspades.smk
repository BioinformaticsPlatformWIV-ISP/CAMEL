from pathlib import Path

from camel.app.camel import Camel
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import gene_detection
from camel.scripts.enterococcuspipeline.snakefile import plasmidspades as plasmidspades_workflow

checkpoint plasmidspades_run:
    """
    Runs plasmid SPAdes on the untrimmed reads.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io',
    output:
        FASTA_Contig = Path(config['working_dir']) / 'plasmidspades' / 'fasta.io',
        INFORMS = Path(config['working_dir']) / plasmidspades_workflow.OUTPUT_PLASMIDSPADES_INFORMS
    params:
        working_dir = Path(config['working_dir']) / 'plasmidspades',
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE',
        spades_options = config.get('assembly', {}).get('spades', {}),
    run:
        from camel.app.tools.spades.spades import SPAdes
        spades = SPAdes(Camel.get_instance())
        fq_dict = SnakePipelineUtils.extracts_fq_input(input.IO, key_pe='FASTQ_PE_1', keys_se=[
            'FASTQ_SE_1', 'FASTQ_SE_2'], key_se='FASTQ_SE_1', drop_empty=True, read_type=params.read_type)
        spades.add_input_files(fq_dict)
        step = Step(rule, spades, Camel.get_instance(), params.working_dir, config)
        spades.update_parameters(plasmid=True)
        spades.update_parameters(**params.spades_options)
        step.run_step()
        if 'FASTA_Contig' in spades.tool_outputs:
            SnakemakeUtils.dump_tool_output(spades, 'FASTA_Contig', output.FASTA_Contig)
        else:
            SnakemakeUtils.dump_object([], output.FASTA_Contig)
        spades.informs['_tag'] = 'plasmidSPAdes'
        SnakemakeUtils.dump_object(spades.informs, output.INFORMS)

rule plasmidspades_quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = lambda wildcards: checkpoints.plasmidspades_run.get().output.FASTA_Contig
    output:
        TSV = Path(config['working_dir']) / 'plasmidspades' / 'quast' / 'tsv.io'
    params:
        running_dir = Path(config['working_dir']) / 'plasmidspades' / 'quast'
    run:
        from camel.app.tools.quast.quast import Quast
        quast = Quast(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(quast, input)
        step = Step(rule, quast, Camel.get_instance(), params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast, output)

rule plasmidspades_quast_informs:
    """
    Extracts the information from the QUAST output file.
    """
    input:
        TSV = rules.plasmidspades_quast.output.TSV
    output:
        INFORMS = Path(config['working_dir']) / 'plasmidspades' / 'quast' / 'informs.io'
    params:
        running_dir = Path(config['working_dir']) / 'plasmidspades' / 'quast'
    run:
        from camel.app.tools.quast.quastinformextractor import QuastInformExtractor
        quast_inform_extractor = QuastInformExtractor(Camel.get_instance())
        SnakemakeUtils.add_pickle_inputs(quast_inform_extractor, input)
        step = Step(rule, quast_inform_extractor, Camel.get_instance(), params.running_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(quast_inform_extractor, output)

rule plasmidspades_assembly_report:
    """
    Runs gene detection with the resfinder database using the plasmidSPAdes assembly.
    """
    input:
        FASTA = lambda wildcards: checkpoints.plasmidspades_run.get().output.FASTA_Contig if plasmidspades_workflow.plasmidspades_successful(checkpoints.plasmidspades_run) else [],
        INFORMS_quast = lambda wildcards: rules.plasmidspades_quast_informs.output.INFORMS if plasmidspades_workflow.plasmidspades_successful(checkpoints.plasmidspades_run) else [],
        INFORMS_spades = lambda wildcards: checkpoints.plasmidspades_run.get().output.INFORMS,
    output:
        HTML = Path(config['working_dir']) / plasmidspades_workflow.OUTPUT_PLASMIDSPADES_REPORT
    params:
        dir_working = Path(config['working_dir']) / 'plasmidspades' / 'report'
    run:
        from camel.app.tools.pipelines.assembly.htmlreporterplasmidspades import HTMLReporterPlasmidSpades
        reporter = HTMLReporterPlasmidSpades(Camel.get_instance())
        input_not_empty = {key: value for key, value in input.items() if len(value) > 0}
        SnakemakeUtils.add_pickle_inputs(reporter, input_not_empty)
        step = Step(rule, reporter, Camel.get_instance(), params.dir_working, config)
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule plasmidspades_gene_detection:
    """
    Performs gene detection with on the plasmidSPAdes assembly.
    """
    input:
        FASTA = lambda wildcards: checkpoints.plasmidspades_run.get().output.FASTA_Contig
    output:
        HTML = Path(config['working_dir']) / 'plasmidspades' / 'gene_detection' / '{db}' / 'html_not_empty.io',
        INFORMS = Path(config['working_dir']) / 'plasmidspades' / 'gene_detection' / '{db}' / 'informs_not_empty.io',
    params:
        dir_working = lambda wildcards: Path(config['working_dir']) / 'plasmidspades' / wildcards.db,
        db_config = lambda wildcards: config['gene_detection'][wildcards.db],
        sample_name = config['sample_name']
    threads: 4
    run:
        import bs4
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.components.workflows.genedetectionwrapper import GeneDetectionWrapper
        from camel.app.components.html.htmlreportsection import HtmlReportSection

        fasta_input = SnakemakeUtils.load_object(input.FASTA)
        wrapper = GeneDetectionWrapper(params.dir_working)
        wrapper.run_workflow_blast(fasta_input[0].path, params.sample_name, params.db_config, int(threads))
        section_updated_title = HtmlReportSection(None)
        soup = bs4.BeautifulSoup(wrapper.output.report_section.to_html(), 'html.parser')

        # Replace header
        soup.find('h3').string = f"{soup.find('h3').text} - plasmidSPades"

        # Fix the path for the included files
        for a in soup.find_all('a'):
            if not a.get('href', '').startswith('gene_detection'):
                continue
            a['href'] = Path('plasmidspades') / a['href']
        for path_file, path_symlink in wrapper.output.report_section.files:
            section_updated_title.add_file(path_file, Path('plasmidspades') / path_symlink)

        # Export HTML
        section_updated_title.add_raw(soup.find('div', class_='report_section').decode_contents())
        SnakemakeUtils.dump_object([ToolIOValue(section_updated_title)], output.HTML)

        # Export informs
        SnakemakeUtils.dump_object([hit.to_table_row() for hit in wrapper.output.detected_hits], output.INFORMS)

rule plasmidspades_gene_detection_select_output:
    """
    Selects the output report for the plasmid spades gene detection.
    - Assembly OK -> Gene detection output
    - Assembly failed -> Empty report
    """
    input:
        HTML = lambda wildcards: rules.plasmidspades_gene_detection.output.HTML if plasmidspades_workflow.plasmidspades_successful(checkpoints.plasmidspades_run) else [],
        INFORMS = lambda wildcards: rules.plasmidspades_gene_detection.output.INFORMS if plasmidspades_workflow.plasmidspades_successful(checkpoints.plasmidspades_run) else [],
        HTML_empty = lambda wildcards: str(Path(config['working_dir'] / gene_detection.OUTPUT_GENE_DETECTION_REPORT_EMPTY)) if not plasmidspades_workflow.plasmidspades_successful(checkpoints.plasmidspades_run) else [],
    output:
        HTML = Path(config['working_dir']) / plasmidspades_workflow.OUTPUT_PLASMIDSPADES_GENE_DETECTION_REPORT,
        INFORMS = Path(config['working_dir']) / 'plasmidspades' / 'gene_detection' / '{db}' / 'informs.io',
    run:
        import shutil
        html_file = input.HTML if len(input.HTML) > 0 else input.HTML_empty
        shutil.copyfile(str(html_file), str(output.HTML))
        informs = SnakemakeUtils.load_object(str(input.INFORMS)) if len(input.INFORMS) > 0 else []
        SnakemakeUtils.dump_object(informs, output.INFORMS)

rule plasmidspades_summary:
    """
    Dumps the summary information from the plasmidSPAdes assembly workflow.
    """
    input:
        INFORMS_quast = lambda wildcards: rules.plasmidspades_quast_informs.output.INFORMS if plasmidspades_workflow.plasmidspades_successful(checkpoints.plasmidspades_run) else [],
        INFORMS_gene_detection = expand(rules.plasmidspades_gene_detection_select_output.output.INFORMS, db=('resfinder', 'ncbi_amr'))
    output:
        TSV = Path(config['working_dir']) / plasmidspades_workflow.OUTPUT_PLASMIDSPADES_SUMMARY
    params:
        running_dir = Path(config['working_dir']) / 'plasmidspades' / 'summary'
    run:
        import json
        quast_informs = SnakemakeUtils.load_object(input.INFORMS_quast) if len(input.INFORMS_quast) > 1 else {}
        summary_data = [
            ('plasmidspades_assembly_n50', quast_informs.get('contig', {}).get('N50')),
            ('plasmidspades_assembly_nb_contigs', quast_informs.get('contig', {}).get('# contigs')),
            ('plasmidspades_assembly_total_length', quast_informs.get('genome', {}).get('Total length'))
        ]
        for path_informs in [Path(x) for x in input.INFORMS_gene_detection]:
            summary_data.append(
                (f'hits_{path_informs.parent.name}', json.dumps(SnakemakeUtils.load_object(str(path_informs)))))

        with open(output.TSV, 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')
