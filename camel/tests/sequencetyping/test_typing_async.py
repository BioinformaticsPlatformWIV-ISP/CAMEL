import unittest

from camel.app.camel import Camel
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.sequence_typing.typeasync import TypeAsync


class TestTypingAsync(CamelTestSuite):
    """
    Tests the asynchronous sequence typing tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('typing')
    input_db_nucl = test_file_dir / 'scheme_mlst_neisseria'
    input_db_protein = test_file_dir / 'scheme_pora_neisseria'
    input_db_mixed = test_file_dir / 'scheme_fhbp_neisseria'
    input_fasta = test_file_dir / 'neisseria_mc58.fasta'
    input_typing_reads = {
        'illumina': [test_file_dir / 'S15BD05018_S58_L001_1.fastq', test_file_dir / 'S15BD05018_S58_L001_2.fastq'],
        'iontorrent': [test_file_dir / 'ERR1447913_ds.fastq'],
        'nanopore': [test_file_dir / 'ERR2259087.fastq.gz']
    }

    def test_typing_async_blast_nucl(self) -> None:
        """
        Tests typing with BLAST and a nucleotide scheme.
        :return: None
        """
        typing_async = TypeAsync(Camel.get_instance())
        typing_async.add_input_files({
            'DIR': [ToolIODirectory(TestTypingAsync.input_db_nucl)],
            'FASTA': [ToolIOFile(TestTypingAsync.input_fasta)],
        })
        typing_async.run(self.running_dir)
        self.assertGreater(len(typing_async.tool_outputs), 0)
        self.assertIn('_command', typing_async.informs)

    def test_typing_async_blast_protein(self) -> None:
        """
        Tests typing with BLAST and a protein scheme.
        :return: None
        """
        typing_async = TypeAsync(Camel.get_instance())
        typing_async.add_input_files({
            'DIR': [ToolIODirectory(TestTypingAsync.input_db_protein)],
            'FASTA': [ToolIOFile(TestTypingAsync.input_fasta)],
        })
        typing_async.run(self.running_dir)
        self.assertGreater(len(typing_async.tool_outputs), 0)
        self.assertIn('_command', typing_async.informs)

    def test_typing_async_blast_mixed(self) -> None:
        """
        Tests typing with BLAST and a mixed scheme.
        :return: None
        """
        typing_async = TypeAsync(Camel.get_instance())
        typing_async.add_input_files({
            'DIR': [ToolIODirectory(TestTypingAsync.input_db_protein)],
            'FASTA': [ToolIOFile(TestTypingAsync.input_fasta)],
        })
        typing_async.run(self.running_dir)
        self.assertGreater(len(typing_async.tool_outputs), 0)
        self.assertIn('_command', typing_async.informs)


if __name__ == '__main__':
    unittest.main()
