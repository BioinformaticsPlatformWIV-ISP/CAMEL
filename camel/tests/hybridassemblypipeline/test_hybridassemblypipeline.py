import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.hybridassemblypipeline.mainhybridassemblypipeline import MainHybridAssemblyPipeline
from camel.tests import longRunningTest


class TestHybridAssemblyPipeline(CamelTestSuite):
    """
    Class to test the hybrid assembly pipeline.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('hybridassembly')
    FASTQ_1 = test_file_dir / 'pilsner_1.fq.gz'
    FASTQ_2 = test_file_dir / 'pilsner_2.fq.gz'
    FASTQ_SE = test_file_dir / 'pilsner_ont.fq.gz'

    @longRunningTest()
    def test_hybrid_assembly(self) -> None:
        """
        Tests the hybrid assembly pipeline.
        :return: None
        """
        path_report_out = self.running_dir / 'output.tsv'
        args = [
            '--fastq-pe', str(TestHybridAssemblyPipeline.FASTQ_1), str(TestHybridAssemblyPipeline.FASTQ_2),
            '--fastq-se', str(TestHybridAssemblyPipeline.FASTQ_SE),
            '--working-dir', str(self.running_dir),
            '--threads', '16',
            '--expected-species', 'bacillus_velezensis',
            '--expected-gc-content', '45',
            '--expected-genome-size', '4.5m',
            '--ont-qual', 'nano-corr'
        ]
        main = MainHybridAssemblyPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
