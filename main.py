import logging

from app.camel import Camel
from app.io.tool_io_file import ToolIOFile
from app.tools.blast.blastn import Blastn
from app.tools.fastqc.fastqc import FastQC
from app.tools.srst2.srst2_mlst import Srst2Mlst
from app.tools.trimmomatic.trimmomatic import Trimmomatic


def run_example_fastqc(camel):
    """
    Runs the FastQC example.

    Demonstrate:
        - Missing input
        - Non existing input
        - Default parameters
        - Updating parameters
        - Non existing parameter
        - Tool run directory
        - Non existing run directory
    :return: None
    """
    fastq_files = [ToolIOFile('/home/bebog/sync/camel_2.0/data/r_1.fq'),
                   ToolIOFile('/home/bebog/sync/camel_2.0/data/r_2.fq')]

    fastqc = FastQC(camel)
    fastqc.update_parameters(non_existing_parameter=16, quiet=False)
    for param_name in fastqc._parameters:
        print('{}: {}'.format(param_name, str(fastqc._parameters[param_name])))

    fastqc.add_input_files({'BAM': fastq_files})
    fastqc.update_parameters(quiet=False, threads=10)
    fastqc.run('fastqc_out')

    # Checking the output
    print(fastqc.stdout)
    print(fastqc.stderr)
    fastqc.display()

    print([f.path for f in fastqc.get_outputs('HTML')])


def run_example_trimmomatic(camel):
    """
    Runs the trimmomatic example.
    :return: None
    """
    fastq_files = [ToolIOFile('/home/bebog/sync/camel_2.0/data/r_1.fq'),
                   ToolIOFile('/home/bebog/sync/camel_2.0/data/r_2.fq')]

    trimmomatic = Trimmomatic(camel)
    trimmomatic.add_input_files({'FASTQ_PE': fastq_files})
    trimmomatic.run('trimmomatic_out')
    trimmomatic.display()
    print([f.path for f in trimmomatic.get_outputs('FASTQ_PE')])


def run_example_srst2(camel):
    """
    Runs the SRST2 example.
    Demonstrate:
        - Logging to main log (camel.log) and tool specific log
    :return: None
    """
    logging.info('Starting the SRST2 example')
    fastq_files = [ToolIOFile('/data/testdata/neisseria/160511_miseq/2012-79_S7_L001_R1_001.fastq.gz'),
                   ToolIOFile('/data/testdata/neisseria/160511_miseq/2012-79_S7_L001_R2_001.fastq.gz')]
    db_sequences = [ToolIOFile('/home/bebog/sync/camel_2.0/data/mlst_sequences.fasta')]

    srst2 = Srst2Mlst(camel)
    srst2.add_input_files({'FASTQ_PE': fastq_files, 'FASTA': db_sequences})
    srst2.update_parameters(delimiter='_', threads=24)
    srst2.run('srst2_out')


def run_example_blast(camel):
    """
    Runs the Blast example.
    :param camel: Camel camel instance.
    :return: None
    """
    blastn = Blastn(camel)
    fasta_file_query = ToolIOFile('/home/bebog/sync/camel_2.0/data/contigs.fasta')
    fasta_file_subject = ToolIOFile('/home/bebog/sync/camel_2.0/data/mlst_sequences.fasta')

    blastn.add_input_files({'FASTA': [fasta_file_query], 'FASTA_Subject': [fasta_file_subject]})
    blastn.update_parameters(output_format=6)
    blastn.run('blastn_out')
    print(blastn._tool_outputs)
    print(blastn.get_outputs('TSV'))

if __name__ == '__main__':
    c = Camel('/home/bebog/config/db.yml')
    run_example_fastqc(c)
    # run_example_trimmomatic(c)
    # run_example_srst2(c)
    # run_example_blast(c)
