from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import assembly
from camel.scripts.salmonellapipeline.snakefile import serotyping_seqsero2

camel = Camel.get_instance()

rule serotyping_seqsero2_run_fasta:
    """
    Runs SeqSero2 to perform serotyping on FASTA input.
    Note: The detection method is always (genome assembly) 'kmer' for FASTA input.
    """
    input:
        FASTA = Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA
    output:
        TXT = Path(config['working_dir']) / 'serotyping_seqsero2' / 'serotyping_seqsero2_kmer' / 'SeqSero2_result.io',
        INFORMS = Path(config['working_dir']) / 'serotyping_seqsero2' / 'serotyping_seqsero2_kmer' / 'informs.io'
    params:
        mode = 'kmer',
        running_dir = Path(config['working_dir']) / 'serotyping_seqsero2' / f'serotyping_seqsero2_kmer',
        db_path_seqsero2 = config['serotyping']['seqsero2']['path']
    run:
        from camel.app.tools.pipelines.salmonella.seqsero2 import SeqSero2

        seqsero_tool = SeqSero2(camel)
        seqsero_tool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path_seqsero2)))]})
        SnakemakeUtils.add_pickle_input(seqsero_tool, 'FASTA', Path(input.FASTA))
        seqsero_tool.update_parameters(mode=str(params.mode))
        step = Step(str(rule), seqsero_tool, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqsero_tool, output)

rule serotyping_seqsero2_run_fastq:
    """
    Runs SeqSero2 to perform serotyping on FASTQ input.
    Note: The detection method can be 'allele' or 'kmerread' for FASTQ input.
    """
    input:
        IO = Path(config['working_dir']) / 'fq_dict.io'
    output:
        TXT = Path(config['working_dir']) / 'serotyping_seqsero2' / 'serotyping_seqsero2_{mode}' / 'SeqSero2_result.io',
        INFORMS = Path(config['working_dir']) / 'serotyping_seqsero2' / 'serotyping_seqsero2_{mode}' / 'informs.io'
    params:
        mode = lambda wildcards: wildcards.mode,
        running_dir = lambda wildcards: Path(config['working_dir']) / 'serotyping_seqsero2' / f'serotyping_seqsero2_{wildcards.mode}',
        db_path_seqsero2 = config['serotyping']['seqsero2']['path']
    run:
        from camel.app.tools.pipelines.salmonella.seqsero2 import SeqSero2

        seqsero_tool = SeqSero2(camel)
        seqsero_tool.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path_seqsero2)))]})
        if config['input_type'] == 'illumina':
            seqsero_tool.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        elif config['input_type'] == 'ont':
            seqsero_tool.add_input_files(SnakePipelineUtils.extracts_fq_input(
                Path(input.IO), key_se='FASTQ_ONT', read_type='SE'))
        seqsero_tool.update_parameters(mode=str(params.mode))
        step = Step(str(rule), seqsero_tool, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(seqsero_tool, output)

rule serotyping_seqsero2_dump_summary_info:
    """
    This rule creates a simple output summary for the SeqSero2 serotyping analysis.
    """
    input:
        TXT_seqsero2_kmer = rules.serotyping_seqsero2_run_fasta.output.TXT,
        INFORMS_seqsero2_kmer = rules.serotyping_seqsero2_run_fasta.output.INFORMS,
        TXT_seqsero2_allele = lambda wildcards: str(rules.serotyping_seqsero2_run_fastq.output.TXT).format(mode='allele') if config['input_type'] == 'illumina' else [],  # exclude fasta, ONT & hybrid
        INFORMS_seqsero2_allele = lambda wildcards: str(rules.serotyping_seqsero2_run_fastq.output.INFORMS).format(mode='allele') if config['input_type'] == 'illumina' else [],  # exclude fasta, ONT & hybrid
        TXT_seqsero2_kmerread = lambda wildcards: str(rules.serotyping_seqsero2_run_fastq.output.TXT).format(mode='kmerread') if config['input_type'] in ('ont', 'illumina') else [],
        INFORMS_seqsero2_kmerread = lambda wildcards: str(rules.serotyping_seqsero2_run_fastq.output.INFORMS).format(mode='kmerread') if config['input_type'] in ('ont', 'illumina') else []
    output:
        VAL_TSV_seqsero2 = Path(config['working_dir']) / 'serotyping_seqsero2' / 'summary_out_seqsero2.tsv'
    params:
        running_dir = Path(config['working_dir']) / 'serotyping_seqsero2'
    run:
        # parse obligate SeqSero2 output
        informs_seqsero2_kmer = SnakemakeUtils.load_object(Path(str(input.INFORMS_seqsero2_kmer)))
        tsv_results = serotyping_seqsero2.seqsero2_output_parser(SnakemakeUtils.load_object(Path(str(input.TXT_seqsero2_kmer)))[0].path, 'seqsero2_kmer')
        with Path(output.VAL_TSV_seqsero2).open('w') as handle:
            handle.writelines(item + '\n' for item in tsv_results)

        # parse facultative SeqSero2 output
        if input.TXT_seqsero2_allele:
            files = [(SnakemakeUtils.load_object(Path(str(input.TXT_seqsero2_allele)))[0].path, 'seqsero2_allele')]
            if input.TXT_seqsero2_kmerread:  # kmerread output can only be present if allele output is present
                files.append((SnakemakeUtils.load_object(Path(str(input.TXT_seqsero2_kmerread)))[0].path, 'seqsero2_kmerread'))
            for args_tuple in files:
                tsv_results = serotyping_seqsero2.seqsero2_output_parser(*args_tuple)
                with Path(output.VAL_TSV_seqsero2).open('a') as handle:
                    for item in tsv_results:
                        handle.write(item + '\n')

        with Path(output.VAL_TSV_seqsero2).open('a') as handle:
            handle.write(f"seqsero2_tool_version\t{informs_seqsero2_kmer['_name']}\n")
            handle.write(f"seqsero2_db_version\t{informs_seqsero2_kmer['last_update_date']}\n")

rule serotyping_seqsero2_report:
    """
    Creates the HTML report for SeqSero2.
    If Illumina FASTQ data is available, all three modes are executed.
    If ONT FASTQ data is available, the allele mode is not executed because although the tool supposedly provides 
    support, it doesn't actually work.
    For FASTA and hybrid input, only the FASTA input is used.
    """
    input:
        TXT_seqsero2_kmer = rules.serotyping_seqsero2_run_fasta.output.TXT,
        TXT_seqsero2_allele = lambda wildcards: str(rules.serotyping_seqsero2_run_fastq.output.TXT).format(mode='allele') if config['input_type'] == 'illumina' else [],  # exclude fasta, ONT & hybrid
        TXT_seqsero2_kmerread = lambda wildcards: str(rules.serotyping_seqsero2_run_fastq.output.TXT).format(mode='kmerread') if config['input_type'] in ('ont', 'illumina') else [],
        INFORMS_serotyping_seqsero2 = rules.serotyping_seqsero2_run_fasta.output.INFORMS
    output:
        VAL_HTML = Path(config['working_dir']) / 'serotyping_seqsero2' / 'html_seqsero2.io'
    params:
        running_dir = Path(config['working_dir']) / 'serotyping_seqsero2',
        db_path_seqsero2 = config['serotyping']['seqsero2']['path']
    run:
        from camel.app.tools.pipelines.salmonella.seqsero2reporter import SeqSero2Reporter

        reporter = SeqSero2Reporter(camel)
        reporter.add_input_files({'DIR_seqsero2': [ToolIODirectory(Path(str(params.db_path_seqsero2)))]})
        SnakemakeUtils.add_pickle_inputs(reporter, input, excluded_keys=['TXT_seqsero2_allele', 'TXT_seqsero2_kmerread'])
        if input.TXT_seqsero2_allele:
            SnakemakeUtils.add_pickle_input(reporter, 'TXT_seqsero2_allele', Path(str(input.TXT_seqsero2_allele)))
        if input.TXT_seqsero2_kmerread:
            SnakemakeUtils.add_pickle_input(reporter, 'TXT_seqsero2_kmerread', Path(str(input.TXT_seqsero2_kmerread)))
        step = Step(str(rule), reporter, camel, Path(str(params.running_dir)))
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(reporter, output)

rule serotyping_seqsero2_report_empty:
    """
    Creates an empty HTML report for the SeqSero2 serotyping analysis.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / 'serotyping_seqsero2' / 'html_seqsero2-empty.io'
    params:
        running_dir = Path(config['working_dir']) / 'serotyping_seqsero2'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.pipelines.salmonella.seqsero2reporter import SeqSero2Reporter
        SnakePipelineUtils.create_empty_report_section(SeqSero2Reporter.TITLE, Path(output.VAL_HTML))
