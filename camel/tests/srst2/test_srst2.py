from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.tools.srst2.srst2gene import Srst2Gene
from camel.scripts.srst2.mainsrst2gene import MainSrst2Gene
from camel.scripts.srst2.mainsrst2mlst import MainSrst2Mlst
from camel.tests import longRunningTest


class TestSRST2(CamelTestSuite):
    """
    Tests the SRST2 base tools.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('srst2')
    input_pe_reads = [test_file_dir / f'ERR178148-ds_1.fastq', test_file_dir / f'ERR178148-ds_2.fastq']
    input_gene_db_path = test_file_dir / 'resfinder_clustered_srst2.fasta'
    input_mlst_db_fasta = test_file_dir / 'mlst_warwick_all_seq.fasta'
    input_mlst_db_profiles = test_file_dir / 'profiles.tsv'

    def test_dependencies(self) -> None:
        """
        Tests if the tool dependencies are available.
        :return: None
        """
        srst2 = Srst2Gene(Camel.get_instance())
        for dependency in srst2.dependencies:
            command = Command(f'module load {dependency};')
            command.run(self.running_dir)
            self.assertEqual(command.returncode, 0, f"Dependency '{dependency}' cannot be loaded")

    def test_pe_gene_detection(self) -> None:
        """
        Tests gene detection using PE reads as input
        :return: None
        """
        output_file = self.running_dir / 'srst2_tab_out.tsv'
        args = [
            '--fastq-pe', str(TestSRST2.input_pe_reads[0]), str(TestSRST2.input_pe_reads[1]),
            '--gene-fasta', str(TestSRST2.input_gene_db_path),
            '--working-dir', str(self.running_dir),
            '--output-tsv', str(output_file),
            '--min-coverage', str(50), '--max-divergence', str(25)
        ]
        srst2_gene = MainSrst2Gene(args)
        srst2_gene.run()
        self.assertGreater(output_file.stat().st_size, 0)

    def test_se_gene_detection(self) -> None:
        """
        Tests gene detection using SE reads as input
        :return: None
        """
        output_file = self.running_dir / 'srst2_tab_out.tsv'
        args = [
            '--fastq-se', str(TestSRST2.input_pe_reads[0]),
            '--gene-fasta', str(TestSRST2.input_gene_db_path),
            '--working-dir', str(self.running_dir),
            '--output-tsv', str(output_file),
            '--min-coverage', str(50), '--max-divergence', str(25)
        ]
        srst2_gene = MainSrst2Gene(args)
        srst2_gene.run()
        self.assertGreater(output_file.stat().st_size, 0)

    @longRunningTest()
    def test_se_typing(self) -> None:
        """
        Tests typing using SE input.
        :return: None
        """
        output_file = self.running_dir / 'srst2_tab_out.tsv'
        args = [
            '--fastq-se', str(TestSRST2.input_pe_reads[0]),
            '--locus-fasta', str(TestSRST2.input_mlst_db_fasta),
            '--profiles-tsv', str(TestSRST2.input_mlst_db_profiles),
            '--working-dir', str(self.running_dir),
            '--output-tsv', str(output_file),
            '--max-mismatch', str(4)
        ]
        srst2_mlst = MainSrst2Mlst(args)
        srst2_mlst.run()
        self.assertGreater(output_file.stat().st_size, 0)

    @longRunningTest()
    def test_pe_typing(self) -> None:
        """
        Tests typing using PE input without profiles
        :return: None
        """
        output_file = self.running_dir / 'srst2_tab_out.tsv'
        args = [
            '--fastq-pe', str(TestSRST2.input_pe_reads[0]), str(TestSRST2.input_pe_reads[1]),
            '--locus-fasta', str(TestSRST2.input_mlst_db_fasta),
            '--profiles-tsv', str(TestSRST2.input_mlst_db_profiles),
            '--working-dir', str(self.running_dir),
            '--output-tsv', str(output_file),
            '--max-mismatch', str(4)
        ]
        srst2_mlst = MainSrst2Mlst(args)
        srst2_mlst.run()
        self.assertGreater(output_file.stat().st_size, 0)
