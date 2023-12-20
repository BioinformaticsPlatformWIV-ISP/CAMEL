import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.genedetection.maingenedetection import MainGeneDetection


class TestGeneDetection(CamelTestSuite):
    """
    Tests the gene detection workflow.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()

    # FASTA
    input_fasta = test_file_dir / 'workflows' / 'NC_002695.1.fasta'
    input_fasta_galaxy = test_file_dir / 'workflows' / 'dataset_12.dat'

    # FASTQ
    input_fastq_by_key = {
        'illumina': [test_file_dir / 'gene_detection' / 'illumina' / 'reads_illumina_1.fastq',
                     test_file_dir / 'gene_detection' / 'illumina' / 'reads_illumina_2.fastq'],
        'iontorrent': [test_file_dir / 'gene_detection' / 'iontorrent' / 'reads_iontorrent.fastq'],
        'nanopore': [test_file_dir / 'gene_detection' / 'nanopore' / 'reads_nanopore.fastq']
    }

    # Special
    input_reads_no_hit = [test_file_dir / 'workflows' / 'ecoli_1.fastq',
                          test_file_dir / 'workflows' / 'ecoli_2.fastq']
    input_reads_raw = [test_file_dir / 'gene_detection' / 'reads-ds_1P.fastq',
                       test_file_dir / 'gene_detection' / 'reads-ds_2P.fastq']
    input_reads_raw_galaxy = [test_file_dir / 'workflows' / 'dataset_fwd_11.dat',
                              test_file_dir / 'workflows' / 'dataset_rev_10.dat']
    input_gene_detection_db = test_file_dir / 'gene_detection' / 'db'

    ###############
    # FASTA input #
    ###############
    def test_gene_detection_blast_fasta(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on FASTA input.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestGeneDetection.input_fasta),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_fasta_galaxy(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on FASTA input in Galaxy format.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestGeneDetection.input_fasta_galaxy),
            '--fasta-name', TestGeneDetection.input_fasta.name,
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir)
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_fasta_spaces(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on FASTA input with spaces in the name.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestGeneDetection.input_fasta_galaxy),
            '--fasta-name', '"my reference genome.fasta"',
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir)
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    ###############
    # FASTQ input #
    ###############
    # BLAST
    def test_gene_detection_blast_illumina(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--assembly-kmers', '33,55'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_illumina_fasta_out(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data and saving of the assembly.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        path_fasta_out = self.running_dir / 'report' / 'assembly.fasta'
        args = [
            '--fastq-pe', str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-fasta', str(path_fasta_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--assembly-kmers', '33,55'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(path_fasta_out.stat().st_size, 0)

    def test_gene_detection_blast_illumina_trim(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--trim-reads',
            '--assembly-kmers', '33,55',
            '--adapter', 'TruSeq2'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_iontorrent(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on IonTorrent data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['iontorrent'][0]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'iontorrent',
            '--assembly-kmers', '33,55'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_iontorrent_trim(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on IonTorrent data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['iontorrent'][0]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'iontorrent',
            '--trim-reads',
            '--assembly-kmers', '33,55'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_nanopore(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Nanopore data, including trimming.
        :return: None
        """
        pass

    # KMA
    def test_gene_detection_kma_illumina(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'illumina',
            '--detection-method', 'kma',
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_kma_illumina_trim(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'illumina',
            '--detection-method', 'kma',
            '--adapter', 'TruSeq3',
            '--trim-reads'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_kma_iontorrent(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on IonTorrent data.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['iontorrent'][0]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'iontorrent',
            '--detection-method', 'kma',
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_kma_iontorrent_trim(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on IonTorrent data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['iontorrent'][0]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'iontorrent',
            '--detection-method', 'kma',
            '--trim-reads'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_kma_nanopore(self) -> None:
        """
        Tests the gene detection workflow with KMA-based detection on Nanopore data.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['nanopore'][0]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'nanopore',
            '--detection-method', 'kma'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_kma_nanopore_trim(self) -> None:
        """
        Tests the gene detection workflow with KMA-based detection on Nanopore data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['nanopore'][0]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'nanopore',
            '--detection-method', 'kma',
            '--trim-reads'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_kma_ont(self) -> None:
        """
        Tests the gene detection workflow with KMA-based detection on Nanopore data, with ont preset.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['nanopore'][0]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'nanopore',
            '--detection-method', 'kma',
            '--kma-ont'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_kma_apm(self) -> None:
        """
        Tests the gene detection workflow with KMA-based detection on Illumina data, with apm preset.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'illumina',
            '--detection-method', 'kma',
            '--kma-apm', 'p'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


    # SRST2
    def test_gene_detection_srst2_illumina(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'illumina',
            '--detection-method', 'srst2'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_srst2_illumina_trim(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'illumina',
            '--detection-method', 'srst2',
            '--trim-reads'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_srst2_iontorrent(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['iontorrent'][0]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'iontorrent',
            '--detection-method', 'srst2',
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_srst2_iontorrent_trim(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['iontorrent'][0]),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'iontorrent',
            '--detection-method', 'srst2',
            '--trim-reads'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_srst2_nanopore(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Nanopore data, including trimming.
        :return: None
        """
        pass

    #################
    # Special cases #
    #################
    def test_gene_detection_srst2_illumina_read_names(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--fastq-pe-names', 'MB3984_S29_L001_R1_001.fastq.gz', 'MB3984_S29_L001_R2_001.fastq.gz',
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--read-type', 'illumina',
            '--detection-method', 'srst2'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    # def test_gene_detection_blast_fasta_input(self) -> None:
    #     """
    #     Tests the gene detection main script using blast.
    #     :return: None
    #     """
    #     path_report_out = self.running_dir / 'report' / 'report.html'
    #     args = [
    #         '--fasta', str(TestGeneDetection.input_fasta),
    #         '--database-dir', str(TestGeneDetection.input_gene_detection_db),
    #         '--output-html', str(path_report_out),
    #         '--output-dir', str(path_report_out.parent),
    #         '--working-dir', str(self.running_dir)
    #     ]
    #     main = MainGeneDetection(args)
    #     main.run()
    #     self.assertGreater(path_report_out.stat().st_size, 0)
    #
    # def test_gene_detection_blast_fasta_input_galaxy(self) -> None:
    #     """
    #     Tests the gene detection main script using blast with an input file in the Galaxy name style..
    #     :return: None
    #     """
    #     path_report_out = self.running_dir / 'report' / 'report.html'
    #     args = [
    #         '--fasta', str(TestGeneDetection.input_fasta_galaxy),
    #         '--fasta-name', TestGeneDetection.input_fasta.name,
    #         '--database-dir', str(TestGeneDetection.input_gene_detection_db),
    #         '--output-html', str(path_report_out),
    #         '--output-dir', str(path_report_out.parent),
    #         '--working-dir', str(self.running_dir)
    #     ]
    #     main = MainGeneDetection(args)
    #     main.run()
    #     self.assertGreater(path_report_out.stat().st_size, 0)
    #
    # def test_gene_detection_blast_fasta_input_spaces(self) -> None:
    #     """
    #     Tests the gene detection main script using blast.
    #     :return: None
    #     """
    #     path_report_out = self.running_dir / 'report' / 'report.html'
    #     args = [
    #         '--fasta', str(TestGeneDetection.input_fasta_galaxy),
    #         '--fasta-name', '"my reference genome.fasta"',
    #         '--database-dir', str(TestGeneDetection.input_gene_detection_db),
    #         '--output-html', str(path_report_out),
    #         '--output-dir', str(path_report_out.parent),
    #         '--working-dir', str(self.running_dir)
    #     ]
    #     main = MainGeneDetection(args)
    #     main.run()
    #     self.assertGreater(path_report_out.stat().st_size, 0)
    #
    # @longRunningTest()
    # def test_gene_detection_blast_fastq_input(self) -> None:
    #     """
    #     Tests the gene detection main script using blast.
    #     :return: None
    #     """
    #     path_report_out = self.running_dir / 'report' / 'report.html'
    #     args = [
    #         '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw_galaxy],
    #         '--database-dir', str(TestGeneDetection.input_gene_detection_db),
    #         '--output-html', str(path_report_out),
    #         '--output-dir', str(path_report_out.parent),
    #         '--working-dir', str(self.running_dir),
    #         '--trim-reads'
    #     ]
    #     main = MainGeneDetection(args)
    #     main.run()
    #     self.assertGreater(path_report_out.stat().st_size, 0)
    #
    # @longRunningTest()
    # def test_gene_detection_blast_fastq_input_no_trim(self) -> None:
    #     """
    #     Tests the gene detection main script using blast.
    #     :return: None
    #     """
    #     path_report_out = self.running_dir / 'report' / 'report.html'
    #     args = [
    #         '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw_galaxy],
    #         '--database-dir', str(TestGeneDetection.input_gene_detection_db),
    #         '--output-html', str(path_report_out),
    #         '--output-dir', str(path_report_out.parent),
    #         '--working-dir', str(self.running_dir)
    #     ]
    #     main = MainGeneDetection(args)
    #     main.run()
    #     self.assertGreater(path_report_out.stat().st_size, 0)
    #
    # @longRunningTest()
    # def test_gene_detection_blast_galaxy_trimmomatic(self) -> None:
    #     """
    #     Tests the gene detection main script using blast with output generated by trimmomatic in Galaxy.
    #     :return: None
    #     """
    #     path_report_out = self.running_dir / 'report' / 'report.html'
    #     args = [
    #         '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw_galaxy],
    #         '--fastq-pe-names',
    #         'Trimmomatic on Neisseria_2.fastq (R1 paired)', 'Trimmomatic on Neisseria_2.fastq (R2 paired)',
    #         '--database-dir', str(TestGeneDetection.input_gene_detection_db),
    #         '--output-html', str(path_report_out),
    #         '--output-dir', str(path_report_out.parent),
    #         '--working-dir', str(self.running_dir)
    #     ]
    #     main = MainGeneDetection(args)
    #     main.run()
    #     self.assertGreater(path_report_out.stat().st_size, 0)
    #
    # @longRunningTest()
    # def test_gene_detection_blast_invalid_filenames(self) -> None:
    #     """
    #     Tests the gene detection workflow with invalid filenames.
    #     :return: None
    #     """
    #     path_report_out = self.running_dir / 'report' / 'report.html'
    #     args = [
    #         '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw_galaxy],
    #         '--fastq-pe-names', 'InvalidFilename 123', 'InvalidReverseReads 555',
    #         '--database-dir', str(TestGeneDetection.input_gene_detection_db),
    #         '--output-html', str(path_report_out),
    #         '--output-dir', str(path_report_out.parent),
    #         '--working-dir', str(self.running_dir)
    #     ]
    #     main = MainGeneDetection(args)
    #     main.run()
    #     self.assertGreater(path_report_out.stat().st_size, 0)
    #
    # @longRunningTest()
    # def test_gene_detection_srst2(self) -> None:
    #     """
    #     Tests the gene detection main script using SRST2.
    #     :return: None
    #     """
    #     path_report_out = self.running_dir / 'report' / 'report.html'
    #     args = [
    #         '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw],
    #         '--database-dir', str(TestGeneDetection.input_gene_detection_db),
    #         '--output-html', str(path_report_out),
    #         '--output-dir', str(path_report_out.parent),
    #         '--working-dir', str(self.running_dir),
    #         '--detection-method', 'srst2',
    #         '--trim-reads'
    #     ]
    #     main = MainGeneDetection(args)
    #     main.run()
    #     self.assertGreater(path_report_out.stat().st_size, 0)
    #
    # def test_gene_detection_srst2_galaxy(self) -> None:
    #     """
    #     Tests the gene detection main script using SRST2 with Galaxy style inputs.
    #     :return: None
    #     """
    #     path_report_out = self.running_dir / 'report' / 'report.html'
    #     args = [
    #         '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw_galaxy],
    #         '--fastq-pe-names', *[x.name for x in TestGeneDetection.input_reads_raw],
    #         '--database-dir', str(TestGeneDetection.input_gene_detection_db),
    #         '--output-html', str(path_report_out),
    #         '--output-dir', str(path_report_out.parent),
    #         '--working-dir', str(self.running_dir),
    #         '--detection-method', 'srst2'
    #     ]
    #     main = MainGeneDetection(args)
    #     main.run()
    #     self.assertGreater(path_report_out.stat().st_size, 0)
    #
    # def test_gene_detection_srst2_no_hits(self) -> None:
    #     """
    #     Tests the gene detection main script using SRST2 where no hits are detected.
    #     :return: None
    #     """
    #     path_report_out = self.running_dir / 'report' / 'report.html'
    #     args = [
    #         '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_no_hit],
    #         '--database-dir', str(TestGeneDetection.input_gene_detection_db),
    #         '--output-html', str(path_report_out),
    #         '--output-dir', str(path_report_out.parent),
    #         '--working-dir', str(self.running_dir),
    #         '--detection-method', 'srst2',
    #         '--trim-reads'
    #     ]
    #     main = MainGeneDetection(args)
    #     main.run()
    #     self.assertGreater(path_report_out.stat().st_size, 0)
    #
    # def test_gene_detection_srst2_galaxy_trimmomatic(self) -> None:
    #     """
    #     Tests the gene detection main script using blast with output generated by trimmomatic in Galaxy.
    #     :return: None
    #     """
    #     path_report_out = self.running_dir / 'report' / 'report.html'
    #     args = [
    #         '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw_galaxy],
    #         '--fastq-pe-names',
    #         'Trimmomatic on Neisseria_2.fastq (R1 paired)', 'Trimmomatic on Neisseria_2.fastq (R2 paired)',
    #         '--database-dir', str(TestGeneDetection.input_gene_detection_db),
    #         '--output-html', str(path_report_out),
    #         '--output-dir', str(path_report_out.parent),
    #         '--working-dir', str(self.running_dir),
    #         '--detection-method', 'srst2',
    #     ]
    #     main = MainGeneDetection(args)
    #     main.run()
    #     self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
