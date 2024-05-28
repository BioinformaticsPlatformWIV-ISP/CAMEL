from pathlib import Path

import pandas as pd

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.pipeline.step import Step

camel = Camel.get_instance()


rule all:
    input:
        CSV = expand(str(Path(config['output_dir']) / 'results' / 'full_report' / '{primer_name}_results.xlsx'),
            primer_name=config['primers'].keys())

checkpoint parts:
    """
    This rule takes an input FASTA file, divides it into parts with Seqkit and saves each part as a separate FASTA file.
    In the input you can give in how many parts it should be divided, so that downstream analysis can run the commands
    in parallel.
    """
    input:
        FASTA = Path(config['input']['fasta'])
    output:
        DIR = directory(Path(config['output_dir']) / 'results' / 'parts')
    params:
        parts_size = config['parts_size'],
        working_dir = Path(config['working_dir'])
    run:
        from camel.app.tools.seqkit.seqkitsplit2 import SeqkitSplit2

        seq_kit_split_2 = SeqkitSplit2(Camel.get_instance())
        seq_kit_split_2.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA))]})
        seq_kit_split_2.update_parameters(by_part=int(config['split_fasta_file']),
            output_dir=output.DIR)
        step = Step(str(rule),seq_kit_split_2,camel,params.working_dir,config)
        step.run_step()

rule per_parts:
    """
    This rule processes the FASTA parts. It will use the self-written function FindMatches to find for each primer
    if there are degenerate nucleotides within the primer. Next, it will match the primer sequence with the sequences
    within the FASTA part, allowing a certain percentage of mismatches between the primer sequence and the FASTA. Also,
    in the input the number of nucleotides at the end of the sequences where no mismatches can happen are taken into 
    account. This returns a dataframe saved in a CSV with columns: id (of fasta sequence), start (position within the 
    sequence where match starts), stop (position within the sequence where match stops), dist (the number of mismatches 
    between the primer sequence and the FASTA sequence) and matched (the sequence from the FASTA sequence that matched).
    """
    input:
        FASTA = Path(config['output_dir']) / 'results' / 'parts' / '{parts_nb}.fasta'
    output:
        CSV = Path(config['output_dir']) / 'results' / 'processing' / '{primer_name}' / 'parts_{parts_nb}.csv'
    params:
        working_dir = Path(config['working_dir']),
        fasta_primer_name = lambda wildcards: str(wildcards.primer_name),
        primer = lambda wildcards: config['primers'][str(wildcards.primer_name)][0],
        end_mismatch = config['end_mismatch'],
        perc_mismatch = config['perc_mismatch']
    run:
        from camel.scripts.screened2.tools.per_chunk_find_matches import FindMatches

        fm = FindMatches(camel)
        fm.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA))]})
        fm.update_parameters(output=str(output.CSV),primer=str(params.primer),
            fasta_primer_name=str(params.fasta_primer_name),
            end_mismatch=int(params.end_mismatch),
            perc_mismatch=float(params.perc_mismatch))
        step = Step(str(rule),fm,camel,params.working_dir,config)
        step.run_step()


def aggregate_input(wildcards):
    """
    This function takes a 'wildcards' argument and returns a list of formatted file paths
    :return: list: A list of formatted file paths.
    """
    return [str(rules.per_parts.output.CSV).format(
        parts_nb=fasta_parts.stem,
        primer_name=wildcards.primer_name)

        for fasta_parts in Path(checkpoints.parts.get(**wildcards).output.DIR).iterdir()
        if not fasta_parts.name.startswith('.')]


rule combine_parts:
    """
    This rule concatenates all resulting dataframes (CSV file) from each FASTA (rule per_parts), removes duplicates and 
    outputs a consolidated CSV file for each primer name.
    """
    input:
        CSV_INPUT = aggregate_input
    output:
        CSV_OUTPUT = Path(config['output_dir']) / 'results' / 'combined' / '{primer_name}.csv'
    run:
        # Initialize an empty list to store DataFrames
        all_df = pd.DataFrame()
        # Iterate through each CSV file in the sorted input list
        for csv in [Path(x) for x in sorted(input.CSV_INPUT)]:
            # Read the CSV file and store it as a DataFrame.
            parts = pd.read_csv(csv)
            # Append the DataFrame to the all_df list.
            all_df = pd.concat([all_df, parts])

        # Remove duplicate records from the concatenated DataFrame.
        records_out = all_df.drop_duplicates()
        # Write the consolidated DataFrame to the specified output CSV file.
        records_out.to_csv(output.CSV_OUTPUT,index=False)

rule csv_xlsx_all:
    """
    This rule processes CSV and XLSX files containing primer probe check data, extracts relevant information, 
    and generates an XLSX report summarizing the assay results for the specified primer. Relevant fields are extracted 
    such as 'Sequence', 'Start', 'End', 'Mismatches', and 'Matched Sequence' from the provided CSV and XLSX files. 
    The data to generate statistical insights is analyzed on primer performance, including mismatch occurrences, 
    sequence alignments, and assay success rates. A comprehensive report is generated presenting detailed findings, 
    including match details, mismatches, alignment positions, and assay outcomes for the specified primer.
    """
    input:
        CSV = rules.combine_parts.output.CSV_OUTPUT
    output:
        XLSX_STAT = Path(config['output_dir']) / 'results' / 'stat_df' / '{primer_name}_results_summary.xlsx',
        XLSX_FP = Path(config['output_dir']) / 'results' / 'full_report' / '{primer_name}_results.xlsx'
    params:
        working_dir = Path(config['working_dir']),
        fasta_primer_name = lambda wildcards: str(wildcards.primer_name),
        primer = lambda wildcards: config['primers'][str(wildcards.primer_name)][0],
        end_mismatch = config['end_mismatch'],
        perc_mismatch = config['perc_mismatch']
    run:
        from camel.scripts.screened2.tools.primerprobecheckreporter import \
            PrimerProbeCheckReporterCSV

        reporter = PrimerProbeCheckReporterCSV(camel)
        reporter.add_input_files({'CSV': [ToolIOFile(Path(input.CSV))]})
        reporter.update_parameters(output_stat=str(output.XLSX_STAT),
            output_full_report=str(output.XLSX_FP),primer=str(params.primer),
            fasta_primer_name=str(params.fasta_primer_name),
            end_mismatch=int(params.end_mismatch),
            perc_mismatch=float(params.perc_mismatch))
        step = Step(str(rule), reporter, camel, params.working_dir, config)
        step.run_step()
