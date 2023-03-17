from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Kleborate(Tool):
    """
    Kleborate is a tool to screen genome assemblies of Klebsiella pneumoniae and the Klebsiella pneumoniae species
    complex (KpSC) for:
    * MLST sequence type
    * species (e.g. K. pneumoniae, K. quasipneumoniae, K. variicola, etc.)
    * ICEKp associated virulence loci: yersiniabactin (ybt), colibactin (clb), salmochelin (iro), hypermucoidy (rmpA)
    * virulence plasmid associated loci: salmochelin (iro), aerobactin (iuc), hypermucoidy (rmpA, rmpA2)
    * antimicrobial resistance determinants: acquired genes, SNPs, gene truncations and intrinsic β-lactamases
    * K (capsule) and O antigen (LPS) serotype prediction, via wzi alleles and Kaptive
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Kleborate', '2.3.2', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("FASTA input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        path_out = self.folder / 'kleborate.txt'
        self._command.command = ' '.join([
            self._tool_command,
            f"-a {self._tool_inputs['FASTA'][0].path}",
            f'-o {path_out}',
            *self._build_options()
        ])
        self._execute_command()
        self._parse_output_file(path_out)
        self._tool_outputs['TSV'] = [ToolIOFile(path_out)]

    def _parse_output_file(self, path_out: Path) -> None:
        """
        Parses the output file and stores the results in the informs.
        :param path_out: Output file path
        :return: None
        """
        with path_out.open() as handle:
            for line in handle.readlines():
                parts = line.strip().split('\t')
                self._informs[parts[0]] = '\t'.join(parts[1:]).strip()

    def _check_command_output(self) -> None:
        """
        Checks if the command executed successfully.
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f'Error executing {self.name}: {self._command.stderr}')
