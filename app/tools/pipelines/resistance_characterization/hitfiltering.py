from app.components.blast_hit.blasthitclustering import BlastHitClustering
from app.components.blast_hit.blasthitfiltering import BlastHitFiltering
from app.components.blast_hit.blasthitparser import BlastHitParser
from app.io.tooliovalue import ToolIOValue
from app.tools.tool import Tool


class HitFiltering(Tool):
    """
    Class that filters BLAST output to report the best matches.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super(HitFiltering, self).__init__('Resistance Characterization: Hit Filtering', '0.1', camel)

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        hits = BlastHitParser.parse(self._tool_inputs['TSV'][0].path)
        hits = BlastHitFiltering.filter_percent_identity(hits, float(self._parameters['min_percent_identity'].value))
        hits = BlastHitFiltering.filter_coverage(hits, float(self._parameters['min_coverage'].value))
        self._tool_outputs['VAL_Hits'] = [ToolIOValue(hit) for hit in HitFiltering.__get_best_hit_per_position(hits)]

    def _check_input(self):
        """
        Checks whether the input is correct.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise StandardError("No 'TSV' input found.")
        super(HitFiltering, self)._check_input()

    @staticmethod
    def __get_best_hit_per_position(hits):
        """
        Returns the best hit for each group of overlapping BLAST hits.
        :param hits: Input hits
        :return: Best matching hits
        """
        groups = BlastHitClustering.cluster_overlapping(hits)
        reported_hits = []
        for group in groups:
            reported_hits.extend(BlastHitFiltering.detect_best_hits(group.hits))
        return reported_hits
