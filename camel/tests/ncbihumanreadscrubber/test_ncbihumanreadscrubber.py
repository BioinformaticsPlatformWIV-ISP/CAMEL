import unittest

from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.ncbihumanreadscrubber.ncbihumanreadscrubber import NcbiHumanReadScrubber
from camel.scripts.ncbihumanreadscrubber.mainncbihumanreadscrubber import MainNcbiHumanReadScrubber


class TestNcbiHumanReadScrubber(CamelTestSuite):
    """
    Tests the HRRT tool.
    _nh: files without human reads
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('ncbi_human_read_scrubber')

    path_fq_se = test_file_dir / 'Myco-DRR041783-ds_1_subset.fastq'
    path_fq_ont = test_file_dir / 'minion_reads.fastq'
    path_fq_ont_nh = test_file_dir / 'minion_no_human.fastq'
    path_fasta = test_file_dir / 'bacteria_seq_hum.fasta'
    path_fasta_nh = test_file_dir / 'bacteria_seq_no_hum.fasta'
    fasta_list = [path_fasta, path_fasta_nh]
    files_list = [path_fq_se, path_fq_ont, path_fq_ont_nh]
    output = [True, True, False]

    def test_scrubber_fq_ont(self) -> None:
        """
        Tests the scrubber on a not gzipped single end FASTQ file.
        :return: None
        """
        for file, out in zip(TestNcbiHumanReadScrubber.files_list, TestNcbiHumanReadScrubber.output):
            scrubber = NcbiHumanReadScrubber(self.camel)
            scrubber.add_input_files({
                'FASTQ_SINGLE_GUNZIP': [ToolIOFile(file)]
            })
            scrubber.update_parameters(interleaved='false', export_human_reads='true', outputfile=self.running_dir / 'test_scrubber_output.fastq', outputfile_removed=self.running_dir / 'test_scrubber_reads_removed.fastq')
            #updated in the ncbi.py as excluded
            scrubber.run(self.running_dir)
            self.verify_output_files(scrubber, 'FASTQ_SCRUBBED', 1)
            self.verify_output_files(scrubber, 'FASTQ_REMOVED', 1) if out else self.verify_output_files(scrubber, 'FASTQ_REMOVED', 0)

            # Check that the input file is larger than the output file
            self.assertGreater(
                FastqUtils.count_reads(file),
                FastqUtils.count_reads(scrubber.tool_outputs['FASTQ_SCRUBBED'][0].path)
            ) if out else self.assertEqual(
                FastqUtils.count_reads(file),
                FastqUtils.count_reads(scrubber.tool_outputs['FASTQ_SCRUBBED'][0].path)
            )

            # Check if the informs were added
            self.assertIn('statistics', scrubber.informs)

    def test_scrubbing_paired(self) -> None:
        """
        Tests the NCBI human read scrubbing standalone pipeline.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        for name in ['', 'nh_']:
            args = [
                '--fastq-pe',
                str(TestNcbiHumanReadScrubber.test_file_dir / f'{name}reads_illumina_1.fastq'),
                str(TestNcbiHumanReadScrubber.test_file_dir / f'{name}reads_illumina_2.fastq'),
                '--output-html', str(path_report_out),
                '--output-dir', str(path_report_out.parent),
                '--working-dir', str(self.running_dir),
                '--output-tsv', "None",
                '--input-type', 'illumina',
                '--threads', '2',
                '--export-removed-reads'
            ]
            main = MainNcbiHumanReadScrubber(args)
            main.run()
            self.assertGreater(path_report_out.stat().st_size, 0)

    def test_scrubbing_fasta(self) -> None:
        """
        Tests the NCBI human read scrubbing standalone pipeline.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        for fasta in TestNcbiHumanReadScrubber.fasta_list:
            args = [
                '--fasta',
                str(fasta),
                '--output-html', str(path_report_out),
                '--output-dir', str(path_report_out.parent),
                '--working-dir', str(self.running_dir),
                '--output-tsv', "None",
                '--input-type', 'fasta',
                '--threads', '2',
                '--export-removed-reads'
            ]
            main = MainNcbiHumanReadScrubber(args)
            main.run()
            self.assertGreater(path_report_out.stat().st_size, 0)

if __name__ == '__main__':
    unittest.main()
