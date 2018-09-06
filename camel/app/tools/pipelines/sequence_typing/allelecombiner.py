import os
from camel.app.tools.tool import Tool

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile


class AlleleCombiner(Tool):
    """
    This tool combines detected alleles from different loci of a locus set.

    Input:
        VAL_Hits: List containing the hits for each locus of a scheme.
    Output:
        TSV: Tabular output file containing the hits for each locus of the scheme and hit statistics.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('Typing: Allele Combiner', '0.1', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        output_file_path = os.path.join(self._folder, self._parameters['output_filename'].value)
        with open(output_file_path, 'w') as handle:
            handle.write('\t'.join(self._tool_inputs['VAL_Hits'][0].value.get_table_column_names()))
            handle.write('\n')
            for hit in [t.value for t in self._tool_inputs['VAL_Hits']]:
                handle.write(hit.to_table_row())
                handle.write('\n')
        self._tool_outputs['TSV'] = [ToolIOFile(output_file_path)]

    def _check_input(self) -> None:
        """
        Checks if the provided tool input is valid.
        :return: None
        """
        if 'VAL_Hits' not in self._tool_inputs:
            raise InvalidInputSpecificationError("A list of sequence typing hit objects ('VAL_Hits') is required.")
        super()._check_input()
