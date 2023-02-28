from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.pipelinecombine.mainpipelinecombine import MainPipelineCombine


class TestPipelineCombine(CamelTestSuite):
    """
    Tests the PipelineCombine tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('pipeline_combine')

    def test_pipeline_combine_staphylococcus(self) -> None:
        """
        Tests the pipeline combine main script.
        :return: None
        """
        path_tsv_out = self.running_dir / 'summary_combined.tsv'
        pipe_combine = MainPipelineCombine([
            *[str(x) for x in TestPipelineCombine.test_file_dir.glob('pipe_staphylococcus_*')],
            '--output', str(path_tsv_out),
            '--exclude', 'cgmlst-*,qc-*,input_files',
        ])
        pipe_combine.run()
        self.assertGreater(path_tsv_out.stat().st_size, 0)

    def test_pipeline_combine_staphylococcus_with_formatted_genes(self) -> None:
        """
        Tests the pipeline combine main script.
        :return: None
        """
        path_tsv_out = self.running_dir / 'summary_combined.tsv'
        pipe_combine = MainPipelineCombine([
            *[str(x) for x in TestPipelineCombine.test_file_dir.glob('pipe_staphylococcus_*')],
            '--output', str(path_tsv_out),
            '--exclude', 'cgmlst-*,qc*,input_files,downsampling*,assembly*',
            '--gene-format', '{hit[1]} ({hit[2]}%)'
        ])
        pipe_combine.run()
        self.assertGreater(path_tsv_out.stat().st_size, 0)

    def test_pipeline_combine_staphylococcus_include_exclude(self) -> None:
        """
        Tests the pipeline combine main script.
        :return: None
        """
        path_tsv_out = self.running_dir / 'summary_combined.tsv'
        pipe_combine = MainPipelineCombine([
            *[str(x) for x in TestPipelineCombine.test_file_dir.glob('pipe_staphylococcus_*')],
            '--output', str(path_tsv_out),
            '--exclude', 'mlst-*,cgmlst-*,qc*,input_files,downsampling*,assembly*,vfdb*,trimming*',
            '--include', 'mlst-ST,cgmlst-SAUR0001'
        ])
        pipe_combine.run()
        self.assertGreater(path_tsv_out.stat().st_size, 0)
