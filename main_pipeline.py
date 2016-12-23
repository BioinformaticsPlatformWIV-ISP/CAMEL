# Test script for running pipelines
from app.camel import Camel
from app.io.tool_io_file import ToolIOFile
from pipeline.pipeline import Pipeline



contigs = ToolIOFile('/data/testdir/bebog/module_tests/mod_rc/blast-s.suis-NC_012926/input/s_suis-BM407-genome.fasta')
mlst_seqs = ToolIOFile('/data/blastdb/nucleotide/ResFinder/latest/resfinder.fasta')
fastq_files = [ToolIOFile('/home/bebog/sync/camel_2.0/data/r_1.fq'),
               ToolIOFile('/home/bebog/sync/camel_2.0/data/r_2.fq')]

camel = Camel('/home/bebog/config/db.yml')
p = Pipeline('data/example.yml', camel)
p.set_initial_input({'FASTA': [contigs], 'DB_BLAST': [mlst_seqs], 'FASTQ_PE': fastq_files})
p.add_options({'Blast_formatter': {'output_format': '0'}, 'FastQC': {'threads': '8', 'quiet': False}})
p.run('/data/testdir/bebog/camel2.0/output')
