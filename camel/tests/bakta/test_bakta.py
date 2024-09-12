import unittest
from pathlib import Path
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bakta.bakta import Bakta


class TestBakta(CamelTestSuite):
    """
    Tests the BAKTA tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('bakta')
    input_fasta = test_file_dir / 'E_coli_K12_subset.fasta'
    input_spp = test_file_dir / 'B_subtilis_168_subset.fasta'
    input_meta_fasta = test_file_dir / 'metagenome_subset.fasta'

    def test_bakta(self) -> None:
        """
        Tests the Bakta tool.
        :return: None
        """
        bakta = Bakta(self.camel)
        bakta.add_input_files({
            'FASTA': [ToolIOFile(Path(TestBakta.input_fasta))]
        })
        bakta.update_parameters(threads=8)
        bakta.run(self.running_dir)
        self.verify_output_files(bakta, 'FAA_FILE', 1)
        self.verify_output_files(bakta, 'GFF3_FILE', 1)
        self.verify_output_files(bakta, 'GBFF_FILE', 1)

    def test_bakta_args(self) -> None:
        """
        Tests the Bakta tool with added arguments.
        :return: None
        """
        bakta = Bakta(self.camel)
        bakta.add_input_files({
            'FASTA': [ToolIOFile(Path(TestBakta.input_spp))]
        })
        bakta.update_parameters(genus='bacillus', species='subtilis', strain='168', threads=8)
        bakta.run(self.running_dir)
        self.verify_output_files(bakta, 'FAA_FILE', 1)
        self.verify_output_files(bakta, 'GFF3_FILE', 1)
        self.verify_output_files(bakta, 'GBFF_FILE', 1)

    def test_bakta_meta(self) -> None:
        """
        Tests the Bakta tool in metagenome mode.
        :return: None
        """
        bakta = Bakta(self.camel)
        bakta.add_input_files({
            'FASTA': [ToolIOFile(Path(TestBakta.input_meta_fasta))]
        })
        bakta.update_parameters(metagenome=None, threads=8)
        bakta.run(self.running_dir)
        self.verify_output_files(bakta, 'FAA_FILE', 1)
        self.verify_output_files(bakta, 'GFF3_FILE', 1)
        self.verify_output_files(bakta, 'GBFF_FILE', 1)


if __name__ == '__main__':
    unittest.main()
