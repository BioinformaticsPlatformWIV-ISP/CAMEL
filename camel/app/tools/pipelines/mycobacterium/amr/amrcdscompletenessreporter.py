from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class AMRCDSCompletenessReporter(Tool):
    """
    Class that checks the CDS completeness the genes associated with AMR.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('AMR CDS completeness', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'VAL_HITS' not in self._tool_inputs:
            raise InvalidInputSpecificationError(f"'VAL_HITS' input is required")
        if 'DB' not in self._input_informs:
            raise InvalidInputSpecificationError(f"'DB' informs input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the tool.
        :return: None
        """
        # Retrieve gene detection hits and group them by locus
        gene_detection_hits = [io.value for io in self._tool_inputs['VAL_HITS']]
        hit_by_locus = {h.locus: h for h in gene_detection_hits}

        # Determine which loci are missing
        mapping = self._input_informs['DB']['mapping']
        missing_loci = []
        self._informs['missing_loci'] = []
        for seq_id in mapping.keys():
            locus = mapping.get(seq_id)
            if locus in hit_by_locus:
                continue
            missing_loci.append((f'<i>{locus}</i>', mapping.get_metadata(seq_id, 'antibiotics')))
            self._informs['missing_loci'].append(locus)

        # Create output report section
        section = HtmlReportSection('CDS completeness')
        nb_targets = len(mapping)
        section.add_paragraph(
            f'{nb_targets - len(missing_loci)}/{nb_targets} coding regions in the resistance regions BED file are covered at >90% with >90% identity (screening with KMA).')

        # Add table with missing loci
        if len(missing_loci) > 0:
            section.add_header('Missing loci', 4)
            section.add_table([l for l in sorted(missing_loci)], ['Locus', 'Antibiotic(s)'], [('class', 'data')])
            section.add_paragraph(
                'These CDSs are (partially) missing, and this must be taken into account when interpreting the AMR results.')

        # Set output
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]
