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
    files_list = [path_fq_se, path_fq_ont, path_fq_ont_nh]
    output = [True, True, False]

    def test_scrubber_fq_se(self) -> None:
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
            # Updated in the ncbi.py as excluded
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
        Tests the NCBI human read scrubbing standalone pipeline with illumina files.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_removed_reads = self.running_dir / 'out' / 'human_read_scrubbing'
        # _nh: files without human reads
        hum_reads =  ['', 'nh_', '', 'nh_']
        extension = ['', '', '.gz', '.gz']

        for hr, ext in zip(hum_reads, extension):
            args = [
                '--fastq-pe',
                str(TestNcbiHumanReadScrubber.test_file_dir / f'{hr}reads_illumina_1.fastq{ext}'),
                str(TestNcbiHumanReadScrubber.test_file_dir / f'{hr}reads_illumina_2.fastq{ext}'),
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
            (self.assertTrue(path_removed_reads.exists()) and self.assertGreater(path_removed_reads.stat().st_size, 0)) if hr == '' else not self.assertTrue(path_removed_reads.exists())

    def test_scrubbing_fasta(self) -> None:
        """
        Tests the NCBI human read scrubbing standalone pipeline with fasta files.
        :return: None
        """
        path_report_html = self.running_dir / 'out' / 'report.html'
        path_removed_reads = self.running_dir / 'out' / 'human_read_scrubbing'
        human_reads = ['', 'no_']

        for hum in human_reads:
            args = [
                '--fasta',
                str(TestNcbiHumanReadScrubber.test_file_dir / f'bacteria_seq_{hum}hum.fasta'),
                '--output-html', str(path_report_html),
                '--output-dir', str(path_report_html.parent),
                '--working-dir', str(self.running_dir),
                '--output-tsv', "None",
                '--input-type', 'fasta',
                '--threads', '2',
                '--export-removed-reads'
            ]
            main = MainNcbiHumanReadScrubber(args)
            main.run()
            (self.assertTrue(path_removed_reads.exists()) and self.assertGreater(path_removed_reads.stat().st_size, 0)) if hum == '' else not self.assertTrue(path_removed_reads.exists())
            self.assertGreater(path_report_html.stat().st_size, 0)

    def test_scrubbing_default (self) -> None:
        for fastq, val in zip(TestNcbiHumanReadScrubber.files_list, TestNcbiHumanReadScrubber.output):
            scrubber_def = NcbiHumanReadScrubber(self.camel)
            scrubber_def.add_input_files({
                'FASTQ_SINGLE_GUNZIP': [ToolIOFile(fastq)]
            })
            scrubber_def.update_parameters(interleaved='false', export_human_reads='false',
                                       outputfile=self.running_dir / 'test_scrubber_output.fastq')
            # Updated in the ncbi.py as excluded
            scrubber_def.run(self.running_dir)
            self.verify_output_files(scrubber_def, 'FASTQ_SCRUBBED', 1)
            self.verify_output_files(scrubber_def, 'FASTQ_REMOVED', 0)

            # Check that the input file is larger than the output file
            self.assertEqual(
                FastqUtils.count_reads(fastq),
                FastqUtils.count_reads(scrubber_def.tool_outputs['FASTQ_SCRUBBED'][0].path)
            ) if not val else self.assertGreater(
                FastqUtils.count_reads(fastq),
                FastqUtils.count_reads(scrubber_def.tool_outputs['FASTQ_SCRUBBED'][0].path))

            # Check if the informs were added
            self.assertIn('statistics', scrubber_def.informs)
if __name__ == '__main__':
    unittest.main()
