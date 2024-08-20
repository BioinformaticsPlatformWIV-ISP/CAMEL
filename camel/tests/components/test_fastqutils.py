import filecmp
import gzip
import unittest
from pathlib import Path

from Bio import SeqIO

from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.testing.cameltestsuite import CamelTestSuite


class TestFastqUtils(CamelTestSuite):
    """
    Tests the Fastq utils module.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('components')

    def test_get_sample_name_miseq_fmt(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(Path('/data/temp/Z4686_S31_L001_R1_001.fastq.gz')), 'Z4686')

    def test_get_sample_name_miseq_fmt_no_s(self) -> None:
        """
        Tests the get sample name function for MiSEQ format without the sample number.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(Path('/data/temp/Z4686_L001_R1_001.fastq.gz')), 'Z4686')

    def test_get_sample_name_miseq_fmt_no_s_lowercase(self) -> None:
        """
        Tests the get sample name function for MiSEQ format without the sample number.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(Path('/data/temp/Z4686_L001_R1_001.fastq.gz')), 'Z4686')

    def test_get_sample_name_miseq_fmt_no_s_no_l(self) -> None:
        """
        Tests the get sample name function for MiSEQ format without the sample number and lane number.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(Path('/data/temp/Z4686_R1_001.fastq.gz')), 'Z4686')

    def test_get_sample_name_simple_fmt(self) -> None:
        """
        Tests the get sample name function for the simple format.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(Path('/data/temp/Z4686_1.fastq.gz')), 'Z4686')

    def test_get_sample_name_simple_fmt_alt(self) -> None:
        """
        Tests the get sample name function for the simple format.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(Path('/data/temp/S18BD02705_R1.fastq.gz')), 'S18BD02705')

    def test_get_sample_name_invalid_fmt(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        with self.assertRaises(ValueError):
            FastqUtils.get_sample_name(Path('/data/temp/invalid_name.fastq.gz'))

    def test_get_sample_name_with_dots(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(Path('/data/temp/my-sample.1.fastq.gz')), 'my-sample')

    def test_get_sample_name_fq_ext(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(Path('/data/temp/my-sample_1.fq')), 'my-sample')

    def test_get_sample_name_fq_ext_gzipped(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(Path('/data/temp/my-sample_1.fq.gz')), 'my-sample')

    def test_get_sample_name_with_p(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(Path('/data/temp/my-sample_1P.fastq')), 'my-sample')

    def test_get_sample_name_with_underscore(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(Path('/data/temp/my_sample_1.fastq')), 'my_sample')

    def test_get_sample_name_se(self) -> None:
        """
        Tests the get sample name function for a single end sample name.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(
            Path('/data/temp/my_sample.fastq'), FastqUtils.PATTERN_FQ_SE), 'my_sample')

    def test_get_sample_name_se_gzipped(self) -> None:
        """
        Tests the get sample name function for a single end sample name.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(
            Path('/data/temp/my_sample.fastq.gz'), FastqUtils.PATTERN_FQ_SE), 'my_sample')

    def test_get_sample_name_se_gzipped_with_dashes(self) -> None:
        """
        Tests the get sample name function for a single end sample name.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(
            Path('/temp/my_sample22-ds.fastq.gz'), FastqUtils.PATTERN_FQ_SE), 'my_sample22-ds')

    def test_get_sample_name_parentheses(self) -> None:
        """
        Tests the get sample name function that contains parentheses.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(Path('/data/temp/UI-586(SRR7648453)_1.fastq.gz')), 'UI-586SRR7648453')

    def test_get_sample_name_miseq(self) -> None:
        """
        Tests the get sample name function that contains parentheses.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(
            Path('/data/temp/S20BD03018_S16_L001_R1_001.fastq.gz')), 'S20BD03018')

    def test_count_reads_interleaved(self) -> None:
        """
        Tests function that returns the number of reads.
        :return: None
        """
        self.assertEqual(FastqUtils.count_reads(TestFastqUtils.test_file_dir / 'interleaved.fq'), 4)

    def test_count_reads_interleaved_gzipped(self) -> None:
        """
        Tests function that returns the number of reads.
        :return: None
        """
        self.assertEqual(FastqUtils.count_reads(TestFastqUtils.test_file_dir / 'interleaved.fq.gz'), 4)

    def test_sort_fastq_by_identifier(self) -> None:
        """
        Test the function that sorts a fastq file by identifier.
        :return: None
        """
        infile = TestFastqUtils.test_file_dir / 'interleaved.fq'
        outfile = self.running_dir / 'sorted.fq'
        compfile = TestFastqUtils.test_file_dir / 'interleaved_sorted.fq'
        FastqUtils.sort_fastq_by_identifier(infile, outfile)
        self.assertTrue(filecmp.cmp(outfile, compfile, shallow=False))

    def test_sort_fastq_by_identifier_zipped(self) -> None:
        """
        Test the function that sorts a fastq file by identifier.
        :return: None
        """
        infile = TestFastqUtils.test_file_dir / 'interleaved.fq.gz'
        outfile = self.running_dir / 'sorted.fq'
        compfile = TestFastqUtils.test_file_dir / 'interleaved_sorted.fq'
        FastqUtils.sort_fastq_by_identifier(infile, outfile)
        self.assertTrue(filecmp.cmp(outfile, compfile, shallow=False))

    def test_sort_fastq_by_identifier_zipped_output(self) -> None:
        """
        Test the function that sorts a fastq file by identifier and outputs
        a gzipped file.
        :return: None
        """
        infile = TestFastqUtils.test_file_dir / 'interleaved.fq'
        outfile = self.running_dir / 'sorted.fq.gz'
        FastqUtils.sort_fastq_by_identifier(infile, outfile, gzip_output=True)
        self.assertTrue(FileSystemHelper.is_gzipped(outfile))

    def test_convert_fastqs_to_interleaved_fastq(self) -> None:
        """
        Test the function that converts two fastq files to one interleaved fastq.
        :return: None
        """
        in_1 = TestFastqUtils.test_file_dir / 'fwd_reads.fq'
        in_2 = TestFastqUtils.test_file_dir / 'rev_reads.fq'
        correct = TestFastqUtils.test_file_dir / 'interleaved.fq'
        outfile = self.running_dir / 'int.fq'
        FastqUtils.convert_fastqs_to_interleaved_fastq(in_1, in_2, outfile)
        self.assertTrue(filecmp.cmp(correct, outfile, shallow=False))

    def test_convert_fastqs_to_interleaved_fastq_zipped_input(self) -> None:
        """
        Test the function that converts two fastq files to one interleaved fastq
        where the input is gzipped.
        :return: None
        """
        in_1 = TestFastqUtils.test_file_dir / 'fwd_reads.fq.gz'
        in_2 = TestFastqUtils.test_file_dir / 'rev_reads.fq.gz'
        correct = TestFastqUtils.test_file_dir / 'interleaved.fq'
        outfile = self.running_dir / 'int.fq'
        FastqUtils.convert_fastqs_to_interleaved_fastq(in_1, in_2, outfile)
        self.assertTrue(filecmp.cmp(correct, outfile, shallow=False))

    def test_convert_fastqs_to_interleaved_fastq_zipped_output(self) -> None:
        """
        Test the function that converts two fastq files to one interleaved fastq
        where the output is gzipped.
        :return: None
        """
        in_1 = TestFastqUtils.test_file_dir / 'fwd_reads.fq'
        in_2 = TestFastqUtils.test_file_dir / 'rev_reads.fq'
        outfile = self.running_dir / 'int.fq'
        FastqUtils.convert_fastqs_to_interleaved_fastq(in_1, in_2, outfile, gzip_output=True)
        self.assertTrue(FileSystemHelper.is_gzipped(outfile))

    def test_split_interleaved_fastq(self) -> None:
        """
        Test the function that splits an interleaved fastq file into two separate files.
        :return: None
        """
        infile = TestFastqUtils.test_file_dir / 'interleaved.fq'
        out_1 = self.running_dir / 'fwd.fq'
        out_2 = self.running_dir / 'rev.fq'
        FastqUtils.split_interleaved_fastq(infile, out_1, out_2)
        self.assertTrue(filecmp.cmp(TestFastqUtils.test_file_dir / 'fwd_reads.fq', out_1, shallow=False))
        self.assertTrue(filecmp.cmp(TestFastqUtils.test_file_dir / 'rev_reads.fq', out_2, shallow=False))

    def test_split_interleaved_fastq_zipped_input(self) -> None:
        """
        Test the function that splits an interleaved fastq file into two separate files
        where the input is gzipped.
        :return: None
        """
        infile = TestFastqUtils.test_file_dir / 'interleaved.fq.gz'
        out_1 = self.running_dir / 'fwd.fq'
        out_2 = self.running_dir / 'rev.fq'
        FastqUtils.split_interleaved_fastq(infile, out_1, out_2)
        self.assertTrue(filecmp.cmp(TestFastqUtils.test_file_dir / 'fwd_reads.fq', out_1, shallow=False))
        self.assertTrue(filecmp.cmp(TestFastqUtils.test_file_dir / 'rev_reads.fq', out_2, shallow=False))

    def test_split_interleaved_fastq_zipped_output(self) -> None:
        """
        Test the function that splits an interleaved fastq file into two separate files
        where the output is gzipped.
        :return: None
        """
        infile = TestFastqUtils.test_file_dir / 'interleaved.fq.gz'
        out_1 = self.running_dir / 'fwd.fq'
        out_2 = self.running_dir / 'rev.fq'
        FastqUtils.split_interleaved_fastq(infile, out_1, out_2, gzip_output=True)
        self.assertTrue(FileSystemHelper.is_gzipped(out_1))
        self.assertTrue(FileSystemHelper.is_gzipped(out_2))

    def test_create_paired_end(self) -> None:
        """
        Test the function that creates an interleaved paired end file and a single end file
        from two separate files.
        :return: None
        """
        in_1 = TestFastqUtils.test_file_dir / 'fq-fwd.fq'
        in_2 = TestFastqUtils.test_file_dir / 'fq-rev.fq'
        pe_out = self.running_dir / 'pe_test.fq'
        se_out = self.running_dir / 'se.fq'
        FastqUtils.create_paired_end(in_1, in_2, pe_out, se_out)
        self.assertTrue(filecmp.cmp(TestFastqUtils.test_file_dir / 'fq-pe.fq', pe_out, shallow=False))
        self.assertTrue(filecmp.cmp(TestFastqUtils.test_file_dir / 'fq-se.fq', se_out, shallow=False))

    def test_create_paired_end_gzip_input(self) -> None:
        """
        Test the function that creates an interleaved paired end file and a single end file
        from two separate files that are gzipped.
        :return: None
        """
        in_1 = TestFastqUtils.test_file_dir / 'fq-fwd.fq.gz'
        in_2 = TestFastqUtils.test_file_dir / 'fq-rev.fq.gz'
        pe_out = self.running_dir / 'pe.fq'
        se_out = self.running_dir / 'se.fq'
        FastqUtils.create_paired_end(in_1, in_2, pe_out, se_out)
        self.assertTrue(filecmp.cmp(TestFastqUtils.test_file_dir / 'fq-pe.fq', pe_out, shallow=False))
        self.assertTrue(filecmp.cmp(TestFastqUtils.test_file_dir / 'fq-se.fq', se_out, shallow=False))

    def test_create_paired_end_gzip_output(self) -> None:
        """
        Test the function that creates an interleaved paired end file and a single end file
        from two separate files. Output is gzipped.
        :return: None
        """
        in_1 = TestFastqUtils.test_file_dir / 'fq-fwd.fq'
        in_2 = TestFastqUtils.test_file_dir / 'fq-rev.fq'
        pe_out = self.running_dir / 'pe.fq.gz'
        se_out = self.running_dir / 'se.fq.gz'
        FastqUtils.create_paired_end(in_1, in_2, pe_out, se_out, gzip_output=True)
        with gzip.open(pe_out) as outhandle, gzip.open(TestFastqUtils.test_file_dir / 'fq-pe.fq.gz') as check_handle:
            out = [x for x in outhandle]
            check = [x for x in check_handle]
        self.assertEqual(out, check)

    def test_get_all_read_names(self) -> None:
        """
        Test the function that retrieves the read names from the given fastq file.
        :return: None
        """
        read_names = FastqUtils.get_all_read_names(TestFastqUtils.test_file_dir / 'fq-fwd.fq')
        self.assertEqual(read_names, {'M04115:7:000000000-AMA6W:1:1101:18190:1854', 'M04115:7:000000000-AMA6W:1:1101:23504:1916',
                                      'M04115:7:000000000-AMA6W:1:1101:23856:1955', 'M04115:7:000000000-AMA6W:1:1101:9914:1907'})

    def test_get_all_read_names_gzipped_input(self) -> None:
        """
        Test the function that retrieves the read names from the given gzipped fastq file.
        :return: None
        """
        read_names = FastqUtils.get_all_read_names(TestFastqUtils.test_file_dir / 'fq-fwd.fq.gz')
        self.assertEqual(read_names, {'M04115:7:000000000-AMA6W:1:1101:18190:1854', 'M04115:7:000000000-AMA6W:1:1101:23504:1916',
                                      'M04115:7:000000000-AMA6W:1:1101:23856:1955', 'M04115:7:000000000-AMA6W:1:1101:9914:1907'})

    def test_process_paired_end(self) -> None:
        """
        Test the function that creates paired and single end files from multiple input files.
        :return: None
        """
        fwd_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-fwd_part_1.fq'
        fwd_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-fwd_part_2.fq'
        rev_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-rev_part_1.fq'
        rev_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-rev_part_2.fq'
        se_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-se_random_1.fq'
        se_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-se_random_2.fq'
        fwd_out = self.running_dir / 'fwd.fq'
        rev_out = self.running_dir / 'rev.fq'
        se_out = self.running_dir / 'se.fq'
        FastqUtils.process_paired_end([fwd_1, fwd_2], [rev_1, rev_2], [se_1, se_2], fwd_out, rev_out, se_out)
        fwd_reads = {(record.id, record.seq) for record in SeqIO.parse(fwd_out, 'fastq')}
        rev_reads = {(record.id, record.seq) for record in SeqIO.parse(rev_out, 'fastq')}
        se_reads = {(record.id, record.seq) for record in SeqIO.parse(se_out, 'fastq')}
        fwd_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fwd_check.fq', 'fastq')}
        rev_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'rev_check.fq', 'fastq')}
        se_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'se_check.fq', 'fastq')}
        self.assertSetEqual(fwd_reads, fwd_reads_check)
        self.assertSetEqual(rev_reads, rev_reads_check)
        self.assertSetEqual(se_reads, se_reads_check)

    def test_process_paired_end_gzip_input(self) -> None:
        """
        Test the function that creates paired and single end files from multiple gzipped input files.
        :return: None
        """
        fwd_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-fwd_part_1.fq.gz'
        fwd_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-fwd_part_2.fq.gz'
        rev_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-rev_part_1.fq.gz'
        rev_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-rev_part_2.fq.gz'
        se_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-se_random_1.fq.gz'
        se_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-se_random_2.fq.gz'
        fwd_out = self.running_dir / 'fwd.fq'
        rev_out = self.running_dir / 'rev.fq'
        se_out = self.running_dir / 'se.fq'
        FastqUtils.process_paired_end([fwd_1, fwd_2], [rev_1, rev_2], [se_1, se_2], fwd_out, rev_out, se_out)
        fwd_reads = {(record.id, record.seq) for record in SeqIO.parse(fwd_out, 'fastq')}
        rev_reads = {(record.id, record.seq) for record in SeqIO.parse(rev_out, 'fastq')}
        se_reads = {(record.id, record.seq) for record in SeqIO.parse(se_out, 'fastq')}
        fwd_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fwd_check.fq', 'fastq')}
        rev_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'rev_check.fq', 'fastq')}
        se_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'se_check.fq', 'fastq')}
        self.assertSetEqual(fwd_reads, fwd_reads_check)
        self.assertSetEqual(rev_reads, rev_reads_check)
        self.assertSetEqual(se_reads, se_reads_check)

    def test_process_paired_end_gzip_output(self) -> None:
        """
        Test the function that creates paired and single end files from multiple input files.
        Output is gzipped
        :return: None
        """
        fwd_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-fwd_part_1.fq'
        fwd_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-fwd_part_2.fq'
        rev_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-rev_part_1.fq'
        rev_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-rev_part_2.fq'
        se_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-se_random_1.fq'
        se_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-se_random_2.fq'
        fwd_out = self.running_dir / 'fwd.fq.gz'
        rev_out = self.running_dir / 'rev.fq.gz'
        se_out = self.running_dir / 'se.fq.gz'
        FastqUtils.process_paired_end([fwd_1, fwd_2], [rev_1, rev_2], [se_1, se_2], fwd_out, rev_out, se_out, gzip_output=True)
        fwd_reads = {(record.id, record.seq) for record in SeqIO.parse(gzip.open(fwd_out, 'rt'), 'fastq')}
        rev_reads = {(record.id, record.seq) for record in SeqIO.parse(gzip.open(rev_out, 'rt'), 'fastq')}
        se_reads = {(record.id, record.seq) for record in SeqIO.parse(gzip.open(se_out, 'rt'), 'fastq')}
        fwd_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fwd_check.fq', 'fastq')}
        rev_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'rev_check.fq', 'fastq')}
        se_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'se_check.fq', 'fastq')}
        self.assertSetEqual(fwd_reads, fwd_reads_check)
        self.assertSetEqual(rev_reads, rev_reads_check)
        self.assertSetEqual(se_reads, se_reads_check)

    def test_process_paired_end_se(self) -> None:
        """
        Test the function that creates paired and single end files from multiple input files.
        :return: None
        """
        fwd_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-fwd_part_1.fq'
        fwd_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-fwd_part_2.fq'
        rev_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-rev_part_1.fq'
        rev_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-rev_part_2.fq'
        se_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-se_random_1.fq'
        se_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-se_random_2.fq'
        fwd_out = self.running_dir / 'fwd.fq'
        rev_out = self.running_dir / 'rev.fq'
        se_1_out = self.running_dir / 'se_1.fq'
        se_2_out = self.running_dir / 'se_2.fq'
        FastqUtils.process_paired_end_se([fwd_1, fwd_2], [rev_1, rev_2], [se_1], [se_2], fwd_out, rev_out, se_1_out, se_2_out)
        fwd_reads = {(record.id, record.seq) for record in SeqIO.parse(fwd_out, 'fastq')}
        rev_reads = {(record.id, record.seq) for record in SeqIO.parse(rev_out, 'fastq')}
        se_1_reads = {(record.id, record.seq) for record in SeqIO.parse(se_1_out, 'fastq')}
        se_2_reads = {(record.id, record.seq) for record in SeqIO.parse(se_2_out, 'fastq')}
        fwd_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fwd_check.fq', 'fastq')}
        rev_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'rev_check.fq', 'fastq')}
        se_1_reads_check = {(record.id, record.seq) for record in SeqIO.parse(se_1, 'fastq')}
        se_2_reads_check = {(record.id, record.seq) for record in SeqIO.parse(se_2, 'fastq')}
        self.assertSetEqual(fwd_reads, fwd_reads_check)
        self.assertSetEqual(rev_reads, rev_reads_check)
        self.assertSetEqual(se_1_reads, se_1_reads_check)
        self.assertSetEqual(se_2_reads, se_2_reads_check)

    def test_process_paired_end_se_gzip_input(self) -> None:
        """
        Test the function that creates paired and single end files from multiple gzipped input files.
        :return: None
        """
        fwd_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-fwd_part_1.fq.gz'
        fwd_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-fwd_part_2.fq.gz'
        rev_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-rev_part_1.fq.gz'
        rev_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-rev_part_2.fq.gz'
        se_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-se_random_1.fq.gz'
        se_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-se_random_2.fq.gz'
        fwd_out = self.running_dir / 'fwd.fq'
        rev_out = self.running_dir / 'rev.fq'
        se_1_out = self.running_dir / 'se_1.fq'
        se_2_out = self.running_dir / 'se_2.fq'
        FastqUtils.process_paired_end_se([fwd_1, fwd_2], [rev_1, rev_2], [se_1], [se_2], fwd_out, rev_out, se_1_out, se_2_out)
        fwd_reads = {(record.id, record.seq) for record in SeqIO.parse(fwd_out, 'fastq')}
        rev_reads = {(record.id, record.seq) for record in SeqIO.parse(rev_out, 'fastq')}
        se_1_reads = {(record.id, record.seq) for record in SeqIO.parse(se_1_out, 'fastq')}
        se_2_reads = {(record.id, record.seq) for record in SeqIO.parse(se_2_out, 'fastq')}
        fwd_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fwd_check.fq', 'fastq')}
        rev_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'rev_check.fq', 'fastq')}
        se_1_reads_check = {(record.id, record.seq) for record in SeqIO.parse(gzip.open(se_1, 'rt'), 'fastq')}
        se_2_reads_check = {(record.id, record.seq) for record in SeqIO.parse(gzip.open(se_2, 'rt'), 'fastq')}
        self.assertSetEqual(fwd_reads, fwd_reads_check)
        self.assertSetEqual(rev_reads, rev_reads_check)
        self.assertSetEqual(se_1_reads, se_1_reads_check)
        self.assertSetEqual(se_2_reads, se_2_reads_check)

    def test_process_paired_end_se_gzip_output(self) -> None:
        """
        Test the function that creates paired and single end files from multiple input files.
        :return: None
        """
        fwd_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-fwd_part_1.fq'
        fwd_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-fwd_part_2.fq'
        rev_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-rev_part_1.fq'
        rev_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-rev_part_2.fq'
        se_1 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-se_random_1.fq'
        se_2 = TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fq-se_random_2.fq'
        fwd_out = self.running_dir / 'fwd.fq'
        rev_out = self.running_dir / 'rev.fq'
        se_1_out = self.running_dir / 'se_1.fq'
        se_2_out = self.running_dir / 'se_2.fq'
        FastqUtils.process_paired_end_se([fwd_1, fwd_2], [rev_1, rev_2], [se_1], [se_2], fwd_out, rev_out, se_1_out, se_2_out, gzip_output=True)
        fwd_reads = {(record.id, record.seq) for record in SeqIO.parse(gzip.open(fwd_out, 'rt'), 'fastq')}
        rev_reads = {(record.id, record.seq) for record in SeqIO.parse(gzip.open(rev_out, 'rt'), 'fastq')}
        se_1_reads = {(record.id, record.seq) for record in SeqIO.parse(gzip.open(se_1_out, 'rt'), 'fastq')}
        se_2_reads = {(record.id, record.seq) for record in SeqIO.parse(gzip.open(se_2_out, 'rt'), 'fastq')}
        fwd_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'fwd_check.fq', 'fastq')}
        rev_reads_check = {(record.id, record.seq) for record in SeqIO.parse(TestFastqUtils.test_file_dir / 'paired_end_processing' / 'rev_check.fq', 'fastq')}
        se_1_reads_check = {(record.id, record.seq) for record in SeqIO.parse(se_1, 'fastq')}
        se_2_reads_check = {(record.id, record.seq) for record in SeqIO.parse(se_2, 'fastq')}
        self.assertSetEqual(fwd_reads, fwd_reads_check)
        self.assertSetEqual(rev_reads, rev_reads_check)
        self.assertSetEqual(se_1_reads, se_1_reads_check)
        self.assertSetEqual(se_2_reads, se_2_reads_check)

    def test_count_bases(self) -> None:
        """
        Test the function to count the number of bases in a fastq file
        :return: None
        """
        count = FastqUtils.count_bases(TestFastqUtils.test_file_dir / 'fq-pe.fq')
        self.assertEqual(count, 1502)

    def test_count_bases_gzip_input(self) -> None:
        """
        Test the function to count the number of bases in a gzipped fastq file
        :return: None
        """
        count = FastqUtils.count_bases(TestFastqUtils.test_file_dir / 'fq-pe.fq.gz')
        self.assertEqual(count, 1502)

    def test_is_fastq_with_fastq(self) -> None:
        """
        Tests the function that checks whether the input file is a FASTQ file using a FASTQ file as input.
        :return: None
        """
        input_file = TestFastqUtils.test_file_dir / 'fq-file.fastq'
        is_fastq = FastqUtils.is_fastq(input_file)
        self.assertTrue(is_fastq)

    def test_is_fastq_with_gzip_fastq(self) -> None:
        """
        Tests the function that checks whether the input file is a FASTQ file using a gzipped FASTQ file as input.
        :return: None
        """
        input_file = TestFastqUtils.test_file_dir / 'fq-file.fastq.gz'
        is_fastq = FastqUtils.is_fastq(input_file)
        self.assertTrue(is_fastq)

    def test_is_fastq_with_fasta(self) -> None:
        """
        Tests the function that checks whether the input file is a FASTQ file using a FASTA file as input.
        :return: None
        """
        input_file = TestFastqUtils.test_file_dir / 'toy.fasta'
        is_fastq = FastqUtils.is_fastq(input_file)
        self.assertFalse(is_fastq)


if __name__ == '__main__':
    unittest.main()
