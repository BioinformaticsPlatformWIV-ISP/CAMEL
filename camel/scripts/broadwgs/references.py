from pathlib import Path

_dir_references = Path("ref_input")

FASTA_GENOME = _dir_references / "fasta_reference_human_value.io"
FASTA_GENOME_FILE = _dir_references / "fasta_reference_human_value_file.io"
DICT_GENOME = _dir_references / "dictionary_genome_human.io"
DBSNP = _dir_references / "dbsnp_vcf.io"
KNOWN_INDELS = _dir_references / "known_indels_vcf.io"
CALLING_INTERVALS = _dir_references / "calling_intervals.io"
CONTAMINATION_SITES_UD = _dir_references / "contamination_sites.io"
COVERAGE_INTERVALS = _dir_references / "coverage_interval_list.io"
EVALUATION_INTERVALS = _dir_references / "evaluation_interval_list.io"
