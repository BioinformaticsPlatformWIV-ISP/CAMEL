import os
<<<<<<< HEAD

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool
=======
from app.tools.tool import Tool

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliofile import ToolIOFile
>>>>>>> origin/bebog-sequence_typing


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

    def _execute_tool(self):
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

    def _check_input(self):
        """
        Checks if the provided tool input is valid.
        :return: None
        """
        if 'VAL_Hits' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Hit input is required")
        super()._check_input()
