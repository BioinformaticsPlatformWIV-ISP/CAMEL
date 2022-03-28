from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, List, Optional

from camel.app.camel import Camel
from camel.app.components.blast.blasthitstatistics import BLASTN_OUTPUT_FORMAT
from camel.app.components.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.blast.blastformatter import BlastFormatter
from camel.app.tools.blast.blastn import Blastn
from camel.app.tools.blast.blastx import Blastx
from camel.app.tools.kma.kma import KMA
from camel.app.tools.kma.kmatypinghitextractor import KMATypingHitExtractor
from camel.app.tools.pipelines.sequence_typing.alignmentextractor import AlignmentExtractor
from camel.app.tools.pipelines.sequence_typing.besthitselector import BestHitSelector
from camel.app.tools.srst2.srst2alleledetector import SRST2AlleleDetector


@dataclass(frozen=True, unsafe_hash=True)
class TypingResultHolder:
    """
    Class to store the result of a sequence typing job.
    """
    hit: ToolIOValue
    informs: Dict[str, Any]


def detect_hit_blast(
        dir_working: Path, fasta_in: Path, dir_scheme: Path, locus_metadata: Dict, blastn_task: str = 'megablast',
        threads_per_job: int = 1) -> TypingResultHolder:
    """
    Performs hit detection with BLAST+.
    :param fasta_in: Input FASTA file
    :param dir_scheme: Base directory for the typing scheme
    :param locus_metadata: Metadata for the locus
    :param dir_working: Working directory
    :param blastn_task: blastn '-task' parameter value
    :param threads_per_job: Threads per BLAST job
    """
    # Create working directory
    dir_working = Path(dir_working)
    dir_working.mkdir(parents=True, exist_ok=True)

    # Initialize BLAST class
    if locus_metadata['type'] == 'DNA':
        blast = Blastn(Camel.get_instance())
        blast.update_parameters(task=blastn_task)
    elif locus_metadata['type'] == 'peptide':
        blast = Blastx(Camel.get_instance())
    else:
        raise ValueError(f"Invalid locus type: {locus_metadata['type']}")

    # Run BLAST+
    db_path = dir_scheme / locus_metadata['fasta_path']
    blast.add_input_files({'DB_BLAST': [ToolIOFile(db_path)], 'FASTA': [ToolIOFile(fasta_in)]})
    blast.update_parameters(threads=threads_per_job)
    blast.run(dir_working)

    # TSV generation
    formatter_tsv = BlastFormatter(Camel.get_instance())
    formatter_tsv.update_parameters(output_format=BLASTN_OUTPUT_FORMAT)
    formatter_tsv.add_input_files({'ASN': blast.tool_outputs['ASN']})
    formatter_tsv.run(dir_working)

    # Best hit selection
    hit_selector = BestHitSelector(Camel.get_instance())
    hit_selector.add_input_files({'TSV': formatter_tsv.tool_outputs['TSV']})
    hit_selector.add_input_informs({'locus': locus_metadata})
    hit_selector.run(dir_working)

    # Text alignment generation
    formatter_text = BlastFormatter(Camel.get_instance())
    formatter_text.update_parameters(output_format='0', num_alignments=1000)
    formatter_text.add_input_files({'ASN': blast.tool_outputs['ASN']})
    formatter_text.run(dir_working)

    # Alignment extraction
    extractor = AlignmentExtractor(Camel.get_instance())
    extractor.add_input_files({
        'TXT': formatter_text.tool_outputs['TXT'], 'VAL_Hits': hit_selector.tool_outputs['VAL_Hit']})
    extractor.run(dir_working)

    # Add the alignment to the hit object
    if len(extractor.tool_outputs['TXT']) > 0:
        best_hit = hit_selector.tool_outputs['VAL_Hit'][0].value
        best_hit.alignment_path = extractor.tool_outputs['TXT'][0].path
    return TypingResultHolder(hit=hit_selector.tool_outputs['VAL_Hit'][0], informs=blast.informs)


def detect_hit_kma(dir_working: Path, fastq_in: Dict[str, List[ToolIOFile]], dir_scheme: Path, locus_metadata: Dict,
                   read_type: str = 'illumina', threads_per_job: int = 1) -> TypingResultHolder:

    # Create working directory
    dir_working = Path(dir_working)
    dir_working.mkdir(parents=True, exist_ok=True)

    # Get the KMA database
    dir_locus = (dir_scheme / locus_metadata['fasta_path']).parent
    try:
        db_path = next((dir_locus / 'kma').glob('*.name'))
    except StopIteration:
        raise FileNotFoundError(f"KMA database for locus '{locus_metadata['name']}' ({dir_locus}) not found")

    # Launch KMA
    kma = KMA(Camel.get_instance())
    kma.add_input_files(fastq_in)
    kma.add_input_files({'DB': [ToolIOValue(str(db_path.parent / db_path.stem))]})
    if read_type == 'nanopore':
        kma.update_parameters(bc_nano=None, basecalls='0.7')
    kma.run(dir_working)

    # Extract the best hit
    kma_extractor = KMATypingHitExtractor(Camel.get_instance())
    kma_extractor.add_input_files({'TSV': kma.tool_outputs['TSV']})
    kma_extractor.add_input_informs({'locus': locus_metadata})
    kma_extractor.run(dir_working)
    return TypingResultHolder(hit=kma_extractor.tool_outputs['VAL_hit'][0], informs=kma.informs)


def detect_hit_srst2(
        dir_working: Path, fastq_in: Dict[str, List[ToolIOFile]], dir_scheme: Path, locus_metadata: Dict,
        srst2_options: Optional[Dict] = None, threads_per_job: int = 1) -> TypingResultHolder:
    """
    Performs hit detection with BLAST+.
    :param dir_working: Input directory
    :param fastq_in: Input FASTQ dictionary
    :param dir_scheme: Base directory for the typing scheme
    :param locus_metadata: Metadata for the locus
    :param srst2_options: SRST2 options
    :param threads_per_job: Threads per BLAST job
    """
    # Create working directory
    dir_working = Path(dir_working)
    dir_working.mkdir(parents=True, exist_ok=True)

    # Add input files
    detector = SRST2AlleleDetector(Camel.get_instance())
    detector.add_input_files(fastq_in)
    db_path = dir_scheme / locus_metadata['fasta_path']
    detector.add_input_files({'FASTA': [ToolIOFile(db_path)]})
    detector.add_input_informs({'locus': locus_metadata})

    # Update parameters
    if (srst2_options is not None) and ('max_unaligned_overlap' in srst2_options):
        detector.update_parameters(max_unaligned_overlap=srst2_options['max_unaligned_overlap'])
    if 'FASTQ_PE' in fastq_in:
        fwd_read_path = fastq_in['FASTQ_PE'][0].path
        fwd_designator, rev_designator = SequenceTypingUtils.determine_read_status(fwd_read_path)
        detector.update_parameters(forward_designator=fwd_designator, reverse_designator=rev_designator)
    detector.update_parameters(threads=threads_per_job)

    # Run tool
    detector.run(dir_working)
    return TypingResultHolder(detector.tool_outputs['VAL_Hit'][0], detector.informs)


detection_by_method = {
    'blast': detect_hit_blast,
    'kma': detect_hit_kma,
    'srst2': detect_hit_srst2,
}
