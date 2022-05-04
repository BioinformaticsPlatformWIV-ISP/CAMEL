import itertools
import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.mlstphylogeny.mainmlstphylogeny import MainMLSTPhylogeny


class TestMLSTPhylogeny(CamelTestSuite):
    """
    Tests for the MLST phylogeny tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('mlst_phylogeny')
    dir_dataset_large_blast = test_file_dir / 'dataset_large_blast'
    dir_dataset_large_kma = test_file_dir / 'dataset_large_kma'
    dir_dataset_large_srst2 = test_file_dir / 'dataset_large_srst2'
    dir_dataset_small_blast = test_file_dir / 'dataset_small_blast'
    dir_dataset_small_blast_html = test_file_dir / 'dataset_html_small_blast'

    def test_mlst_phylogeny_small_dataset(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and a small dataset.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        input_tsv = list(itertools.chain.from_iterable(
            [['--input-tsv', str(tsv), tsv.name] for tsv in TestMLSTPhylogeny.dir_dataset_small_blast.iterdir()]))
        args = input_tsv + [
            '--output-html', str(output_html),
            '--output-dir', str(output_html.parent),
            '--dir-working', str(self.running_dir)
        ]
        mlst_tree = MainMLSTPhylogeny(args)
        mlst_tree.run()
        self.assertGreater(output_html.stat().st_size, 0)

    def test_mlst_phylogeny_small_dataset_extra_outputs(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and a small dataset and extra output files.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_tsv_alleles = self.running_dir / 'out' / 'allele_matrix.tsv'
        output_tsv_dist = self.running_dir / 'out' / 'dist_matrix.tsv'
        output_nwk = self.running_dir / 'out' / 'tree.nwk'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        input_tsv = list(itertools.chain.from_iterable(
            [['--input-tsv', str(tsv), tsv.name] for tsv in TestMLSTPhylogeny.dir_dataset_small_blast.iterdir()]))
        args = input_tsv + [
            '--output-html', str(output_html),
            '--output-dir', str(output_html.parent),
            '--output-allele-matrix', str(output_tsv_alleles),
            '--output-dist-matrix', str(output_tsv_dist),
            '--output-tree', str(output_nwk),
            '--dir-working', str(self.running_dir)
        ]
        mlst_tree = MainMLSTPhylogeny(args)
        mlst_tree.run()
        for file in [output_html, output_tsv_alleles, output_tsv_dist, output_tsv_alleles]:
            self.assertGreater(file.stat().st_size, 0)

    def test_mlst_phylogeny_small_dataset_html(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and a small dataset and HTML input.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        input_tsv = list(itertools.chain.from_iterable(
            [['--input-html', str(dir_html)] for dir_html in TestMLSTPhylogeny.dir_dataset_small_blast_html.iterdir()]))
        args = input_tsv + [
            '--output-html', str(output_html),
            '--output-dir', str(output_html.parent),
            '--dir-working', str(self.running_dir)
        ]
        mlst_tree = MainMLSTPhylogeny(args)
        mlst_tree.run()
        self.assertGreater(output_html.stat().st_size, 0)

    def test_mlst_phylogeny_large_dataset(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and a large dataset.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        input_tsv = list(itertools.chain.from_iterable(
            [['--input-tsv', str(tsv), tsv.name] for tsv in TestMLSTPhylogeny.dir_dataset_large_blast.iterdir()]))
        args = input_tsv + [
            '--output-html', str(output_html),
            '--output-dir', str(output_html.parent),
            '--dir-working', str(self.running_dir)
        ]
        mlst_tree = MainMLSTPhylogeny(args)
        mlst_tree.run()
        self.assertGreater(output_html.stat().st_size, 0)

    def test_mlst_phylogeny_large_dataset_kma(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and a large dataset and KMA-based detection.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        input_tsv = list(itertools.chain.from_iterable(
            [['--input-tsv', str(tsv), tsv.name] for tsv in TestMLSTPhylogeny.dir_dataset_large_kma.iterdir()]))
        args = input_tsv + [
            '--output-html', str(output_html),
            '--output-dir', str(output_html.parent),
            '--dir-working', str(self.running_dir),
            '--detection-method', 'kma'
        ]
        mlst_tree = MainMLSTPhylogeny(args)
        mlst_tree.run()
        self.assertGreater(output_html.stat().st_size, 0)

    def test_mlst_phylogeny_large_dataset_srst2(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and a large dataset and SRST2-based detection.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        input_tsv = list(itertools.chain.from_iterable(
            [['--input-tsv', str(tsv), tsv.name] for tsv in TestMLSTPhylogeny.dir_dataset_large_srst2.iterdir()]))
        args = input_tsv + [
            '--output-html', str(output_html),
            '--output-dir', str(output_html.parent),
            '--dir-working', str(self.running_dir),
            '--detection-method', 'srst2'
        ]
        mlst_tree = MainMLSTPhylogeny(args)
        mlst_tree.run()
        self.assertGreater(output_html.stat().st_size, 0)

    def test_mlst_phylogeny_identical_profiles(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and identical profiles for all datasets.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        input_tsv = list(sorted(TestMLSTPhylogeny.dir_dataset_small_blast.iterdir()))[0]
        args = [
            '--input-tsv', str(input_tsv), 'dataset_1.tsv',
            '--input-tsv', str(input_tsv), 'dataset_2.tsv',
            '--input-tsv', str(input_tsv), 'dataset_3.tsv',
            '--input-tsv', str(input_tsv), 'dataset_4.tsv',
            '--output-html', str(output_html),
            '--output-dir', str(output_html.parent),
            '--dir-working', str(self.running_dir)
        ]
        mlst_tree = MainMLSTPhylogeny(args)
        mlst_tree.run()
        self.assertGreater(output_html.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
