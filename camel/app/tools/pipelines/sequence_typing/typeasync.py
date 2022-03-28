import json
import logging
import time
from multiprocessing.pool import ThreadPool
from typing import Dict, Any

from camel.app.camel import Camel
from camel.app.components.sequencetyping import typingasyncutils
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
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
        if (self._parameters['detection_method'] == 'blast') and ('FASTA' not in self._tool_inputs):
            raise InvalidInputSpecificationError('FASTA input is required')
        elif (self._parameters['detection_method'] in ('KMA', 'SRST2')) and ('FASTQ_PE' not in self._tool_inputs):
            raise InvalidInputSpecificationError('FASTQ input is required')
        if 'DIR' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Typing directory input is required ('DIR')")
        super()._check_input()

    def get_typing_parameters(self, locus_info) -> Dict[str, Any]:
        """
        Returns the parameters to pass to the typing jobs.
        """
        # General options
        dict_params = {
            'dir_working': self.folder / locus_info['name_sanitized'],
            'dir_scheme': self._tool_inputs['DIR'][0].path,
            'locus_metadata': locus_info,
            'threads_per_job': 1
        }
        # Add tool specific options
        if self._parameters['detection_method'].value == 'blast':
            dict_params['fasta_in'] = self._tool_inputs['FASTA'][0].path
            dict_params['blastn_task'] = self._parameters['blastn_task'].value
        elif self._parameters['detection_method'].value == 'srst2':
            dict_params['fastq_in'] = {k: v for k, v in self._tool_inputs.items() if k.startswith('FASTQ')}
        elif self._parameters['detection_method'].value == 'kma':
            dict_params['fastq_in'] = {k: v for k, v in self._tool_inputs.items() if k.startswith('FASTQ')}
        return dict_params

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        """
        # Parse scheme metadata
        with (self._tool_inputs['DIR'][0].path / 'scheme_metadata.txt').open() as handle:
            metadata = json.load(handle)
        logging.info(f"{len(metadata['loci'])} loci found")

        # Determine which loci need to be typed
        loci_target = [
            info['name_sanitized'] for info in metadata['loci'] if info['type'] == self._parameters['locus_type'].value]

        # Create jobs
        tp = ThreadPool(int(self._parameters['threads'].value))
        jobs = {}
        for locus_info in metadata['loci']:
            if locus_info['name_sanitized'] not in loci_target:
                continue
            locus_fasta = self._tool_inputs['DIR'][0].path / locus_info['fasta_path']
            jobs[locus_info['name_sanitized']] = tp.apply_async(
                typingasyncutils.detection_by_method[self._parameters['detection_method'].value], (),
                self.get_typing_parameters(locus_info))

        # Type all loci
        result_by_locus_name = {}
        while len(jobs) > 0:
            to_remove = []
            for locus_name, result in jobs.items():
                if result.ready():
                    to_remove.append(locus_name)
                    if not result.successful():
                        logging.error(f'Error while typing: {locus_name}')
                        result.get()
            for locus_name in to_remove:
                result_by_locus_name[locus_name] = jobs.pop(locus_name).get()
            logging.debug(f"Job status: {len(jobs)} left (total: {len(loci_target)})")
            time.sleep(0.5)

        # Set output
        self._tool_outputs['VAL_hits'] = [result_by_locus_name[locus].hit for locus in loci_target]

        # Collect informs (for the first locus)
        for key, value in result_by_locus_name[loci_target[0]].informs.items():
            self._informs[key] = value
