import json
import logging
import time
from multiprocessing.pool import ThreadPool
from pathlib import Path
from typing import Dict

from camel.app.camel import Camel
from camel.app.components.blast.blasthitstatistics import BLASTN_OUTPUT_FORMAT
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.blast.blastformatter import BlastFormatter
from camel.app.tools.blast.blastn import Blastn
from camel.app.tools.blast.blastx import Blastx
from camel.app.tools.pipelines.sequence_typing.alignmentextractor import AlignmentExtractor
from camel.app.tools.pipelines.sequence_typing.besthitselector import BestHitSelector
from camel.app.tools.tool import Tool


class TypeAsync(Tool):
    """
    Performs BLAST-based sequence typing asynchronously.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance.
        """
        super().__init__('Typing async', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        if 'DIR' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Typing directory input is required ('DIR')")
        super()._check_input()

    @staticmethod
    def detect_hit_blast(
            dir_working: Path, fasta_in: Path, dir_scheme: Path, locus_metadata: Dict, blastn_task: str = 'megablast',
            threads_per_job: int = 1) -> ToolIOValue:
        """
        Performs hit detection with BLAST+.
        :param fasta_in: Input FASTA file
        :param dir_scheme: Base directory for the typing scheme
        :param locus_metadata: Metadata for the locus
        :param dir_working: Working directory
        :param blastn_task: blastn '-task' parameter value
        :param threads_per_job: Threads per BLAST job
        """
        import pprint
        pprint.pprint(locus_metadata)

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

        # Collect input files
        db_path = dir_scheme / locus_metadata['fasta_path']

        # Run BLAST+
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
        import pprint
        pprint.pprint(hit_selector.tool_outputs['VAL_Hit'][0].__dict__)
        extractor = AlignmentExtractor(Camel.get_instance())
        extractor.add_input_files({
            'TXT': formatter_text.tool_outputs['TXT'],
            'VAL_Hits': hit_selector.tool_outputs['VAL_Hit']})
        extractor.run(dir_working)

        # Add the alignment to the hit object
        if len(extractor.tool_outputs['TXT']) > 0:
            best_hit = hit_selector.tool_outputs['VAL_Hit'][0].value
            best_hit.alignment_path = extractor.tool_outputs['TXT'][0].path
        return hit_selector.tool_outputs['VAL_Hit'][0]

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        """
        # Parse scheme metadata
        with (self._tool_inputs['DIR'][0].path / 'scheme_metadata.txt').open() as handle:
            metadata = json.load(handle)
        logging.info(f"{len(metadata['loci'])} loci found")

        # Create pool with jobs
        tp = ThreadPool(int(self._parameters['threads'].value))
        jobs = {}
        for locus_info in metadata['loci']:
            locus_fasta = self._tool_inputs['DIR'][0].path / locus_info['fasta_path']
            jobs[locus_info['name_sanitized']] = tp.apply_async(
                TypeAsync.detect_hit_blast, (), {
                    'dir_working': self.folder / locus_info['name_sanitized'],
                    'fasta_in': self._tool_inputs['FASTA'][0].path,
                    'dir_scheme': self._tool_inputs['DIR'][0].path,
                    'locus_metadata': locus_info,
                    'blastn_task': 'megablast',
                    'threads_per_job': 1
                })

        # Type all loci
        hit_by_locus_name = {}
        while len(jobs) > 0:
            to_remove = []
            for locus_name, result in jobs.items():
                if result.ready():
                    to_remove.append(locus_name)
                    if not result.successful():
                        logging.error(f'Error while typing: {locus_name}')
                        result.get()
            for locus_name in to_remove:
                hit_by_locus_name[locus_name] = jobs.pop(locus_name).get()
            logging.debug(f"Job STATUS: {len(jobs)}/{len(metadata['loci'])}")
            time.sleep(0.5)

        # Set output
        self._tool_outputs['VAL_hits'] = [
            hit_by_locus_name[locus_info['name_sanitized']] for locus_info in metadata['loci']]
