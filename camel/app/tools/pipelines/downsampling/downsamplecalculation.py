import json
import statistics
from typing import Any, Dict, List

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class DownsampleCalculation(Tool):
    """
    Calculates the factor for downsampling based on the estimated coverage.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Downsample calculation', '0.1', camel)

    def _check_input(self) -> None:
        """
        Check if the provided tool input is valid.
        :return: None
        """
        if 'stats' not in self._input_informs:
            raise InvalidInputSpecificationError("Stats input is required")
        super()._check_input()

    def __calculate_stats(self, fq_stats: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Calculates the downsampling statistics.
        :param fq_stats: FASTQ stats.
        :return: statistics as a dictionary
        """
        ref_genome_size = int(self._parameters['size_ref_genome'].value)
        cov_est = sum(fq['nb_of_bases'] for fq in fq_stats) / ref_genome_size
        cov_target = float(self._parameters['cov_target'].value)
        downsample_factor = float(f"{cov_target / cov_est:.6f}")
        key_nb_reads = 'nb_read_pairs' if 'is_paired' in self._parameters else 'nb_reads'
        data_out = {
            'total_bases': sum([fq['nb_of_bases'] for fq in fq_stats]),
            'mean_read_length': statistics.mean([
                fq['nb_of_bases'] / fq['nb_of_sequences'] for fq in fq_stats]),
            'coverage_estimated': cov_est,
            'coverage_target': cov_target,
            'downsample_factor': cov_target / cov_est if downsample_factor < 1 else None,
            'size_ref_genome': ref_genome_size,
            f'{key_nb_reads}_in': next(iter(fq_stats))['nb_of_sequences']
        }
        return data_out

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        data_stats = self.__calculate_stats(self._input_informs['stats']['stats'])
        self._informs['stats'] = data_stats

        # Create JSON output file
        path_json_out = self.folder / 'downsample_stats.json'
        with path_json_out.open('w') as handle:
            json.dump(data_stats, handle, indent=2)
        self._tool_outputs['JSON'] = [ToolIOFile(path_json_out)]
