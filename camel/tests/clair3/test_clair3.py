import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.vcf.vcfutils import VCFUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.clair3.clair3 import Clair3
from camel.scripts.variantcalling.clair3.maincalling import MainCalling


class TestClair3(CamelTestSuite):
    """
    Initializes this testing tool
    """

    test_file_dir = Path('/testdata/camel/clair3/')
    FILE_FASTA = ToolIOFile(test_file_dir / 'bsubtilis.fa')
    FILE_BAM_ILLUMINA = ToolIOFile(test_file_dir / 'bsubtilis_illumina.bam')
    FILE_BAM_ONT = ToolIOFile(test_file_dir / 'bsubtilis_ont.bam')

    def test_clair3_illumina(self) -> None:
        """
        actually testing Clair3 on illumina sequencing data
        """
        output_file_vcf = self.running_dir / 'merge_output.vcf.gz'
        clair3 = Clair3(self.camel)
        clair3.add_input_files({'FASTA': [TestClair3.FILE_FASTA], 'BAM': [TestClair3.FILE_BAM_ILLUMINA]})
        clair3.update_parameters(output_path=self.running_dir, haploid_precise=True, no_phasing=True, include_ctgs=True)
        clair3.run(self.running_dir)
        self.verify_output_files(clair3, 'VCF')
        self.assertGreater(VCFUtils.count_variants(clair3.tool_outputs['VCF'][0].path), 0)

    def test_clair3_ont(self) -> None:
        """
        actually testing Clair3 on ONT sequencing data
        """
        output_file_vcf = self.running_dir / 'merge_output.vcf.gz'
        clair3 = Clair3(self.camel)
        clair3.add_input_files({'FASTA': [TestClair3.FILE_FASTA], 'BAM': [TestClair3.FILE_BAM_ONT]})
        clair3.update_parameters(output_path=self.running_dir, platform='ont', haploid_precise=True, no_phasing=True,
                                 include_ctgs=True,
                                 model_path='/usr/local/bin/lmod/clair3/0.1.12/bin/models/ont/')
        clair3.run(self.running_dir)
        self.verify_output_files(clair3, 'VCF')
        self.assertGreater(VCFUtils.count_variants(clair3.tool_outputs['VCF'][0].path), 0)

    def test_variant_calling_ont(self) -> None:
        """
        Tests the variant calling main script.
        :return: None
        """
        output_file_vcf = self.running_dir / 'clair3_output.vcf'
        args = [
            '--bam', str(TestClair3.FILE_BAM_ONT),
            '--reference', str(TestClair3.FILE_FASTA),
            '--working-dir', str(self.running_dir),
            '--output', str(output_file_vcf),
            '--model-path', '/usr/local/bin/lmod/clair3/0.1.12/bin/models/ont',
            '--haploid-precise',
            '--no-phasing',
            '--include-ctgs',
            '--platform', 'ont'
        ]
        main_calling = MainCalling(args)
        main_calling.run()
        self.assertGreater(output_file_vcf.stat().st_size, 0)
        self.assertGreater(VCFUtils.count_variants(output_file_vcf), 0)

    def test_variant_calling_illumina(self) -> None:
        """
        Tests the variant calling main script.
        :return: None
        """
        output_file_vcf = self.running_dir / 'clair3_output.vcf'
        args = [
            '--bam', str(TestClair3.FILE_BAM_ILLUMINA),
            '--reference', str(TestClair3.FILE_FASTA),
            '--working-dir', str(self.running_dir),
            '--output', str(output_file_vcf),
            '--model-path', '/usr/local/bin/lmod/clair3/0.1.12/bin/models/ilmn',
            '--haploid-precise',
            '--no-phasing',
            '--include-ctgs',
            '--platform', 'ilmn'
        ]
        main_calling = MainCalling(args)
        main_calling.run()
        self.assertGreater(output_file_vcf.stat().st_size, 0)
        self.assertGreater(VCFUtils.count_variants(output_file_vcf), 0)


if __name__ == '__main__':
    unittest.main()
